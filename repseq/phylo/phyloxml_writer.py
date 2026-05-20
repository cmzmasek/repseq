"""PhyloXML writer with rich per-leaf annotation.

The previous writer was a one-liner: ``Bio.Phylo.write([tree], path,
"phyloxml")`` after restoring leaf names. The output had a ``<name>``
per leaf and nothing else — no taxonomy, no host/date/country, no
sequence accession, no tool provenance. For a bench scientist trying
to interpret a tree of 200 viral isolates, that means flipping back
and forth to a spreadsheet just to find out which leaf is which.

This module replaces that with a writer that emits, per leaf:

* A formatted ``<name>`` driven by the user-configurable label template
  (``repseq/phylo/labels.py``).
* A ``<taxonomy>`` block with NCBI taxon id + scientific name.
* A ``<sequence>`` block with the GenBank accession (``source="ncbi"``)
  and the original FASTA header as ``<name>``.
* Repseq-namespaced ``<property>`` elements for host, collection_date,
  country, strain, isolate_id, year, and the full 9-rank taxonomy
  ladder (species, subgenus, genus, subfamily, family, suborder,
  order, subclass, class). Empty fields are omitted (no zombie
  ``<property/>`` elements).

And, at tree level:

* A ``<phylogeny>`` ``<name>`` and ``<description>`` capturing run
  prefix, alphabet, MSA tool + version, tree tool + version + model +
  bootstrap settings, and a UTC timestamp.

PhyloXML schema-ordering constraint: per ``<clade>``, child elements
must appear in this order: ``name → branch_length → confidence →
taxonomy → sequence → property → child clades``. We build elements in
that order explicitly rather than rely on insertion order, because
some viewers (Archaeopteryx in particular) fall back to defaults if
the schema order is wrong.

Confidence normalisation: tree-builder support values are scaled to
0-100 integers per the rules below (see ``_normalize_confidence``),
and the ``<confidence>`` element gets a ``type=`` attribute naming the
support metric (``sh_like``, ``ufboot``, ``sh_alrt``, ``bootstrap``).
"""

from __future__ import annotations

import logging
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from Bio import Phylo

from .. import __version__ as REPSEQ_VERSION
from ..models import Sequence
from .labels import (
    _parse_year,
    format_leaf_label,
    labeling_options,
    pick_format_string,
)
from .lca import phyloxml_rank

logger = logging.getLogger(__name__)


# XML namespaces. The phyloXML schema URI is fixed; the `repseq:`
# namespace is just a marker on our custom <property> elements so
# downstream tools can tell our annotations apart from anyone else's.
PHYLOXML_NS = "http://www.phyloxml.org"
XSI_NS = "http://www.w3.org/2001/XMLSchema-instance"
REPSEQ_NS = "https://github.com/cmzmasek/repseq"

# Taxonomy ranks exported per leaf, kept in sync with the 9-rank
# ``_TAX_RANKS`` ladder used by the FASTA bracket tags and the rep
# TSVs (``output/report.py``). Sub-ranks (subgenus, suborder,
# subclass) come only via ``TaxonomyInfo.lineage`` and are commonly
# blank for viruses — empty ones are skipped, so a leaf only shows the
# ranks that actually resolved.
_TAXONOMY_RANKS: tuple[str, ...] = (
    "species",
    "subgenus",
    "genus",
    "subfamily",
    "family",
    "suborder",
    "order",
    "subclass",
    "class",
)

# Per-leaf properties exported under the repseq: namespace.
# Each entry is (property_ref, lookup_key). Lookups resolve against the
# Sequence (isolate metadata) or its TaxonomyInfo (rank names). Empty
# results are skipped so the XML carries no zombie ``<property/>``.
_LEAF_PROPERTIES: list[tuple[str, str]] = [
    ("repseq:host", "host"),
    ("repseq:collection_date", "collection_date"),
    ("repseq:country", "country"),
    ("repseq:strain", "strain"),
    ("repseq:isolate_id", "isolate_id"),
    ("repseq:year", "year"),         # derived from collection_date
] + [(f"repseq:{rank}", rank) for rank in _TAXONOMY_RANKS]  # from taxonomy


# ---------------------------------------------------------------------------
# Field lookups
# ---------------------------------------------------------------------------

def _leaf_property_value(seq: Sequence, key: str) -> Optional[str]:
    """Resolve the value for one ``<property>`` entry on a leaf.

    Returns None for missing/empty values; the caller skips them so
    the XML doesn't carry empty stubs.
    """
    if key == "year":
        return _parse_year(seq.collection_date)
    if key in ("host", "collection_date", "country", "strain", "isolate_id"):
        v = getattr(seq, key, None)
        return str(v) if v else None
    if key in _TAXONOMY_RANKS:
        if seq.taxonomy is None:
            return None
        v = seq.taxonomy.get_rank(key)
        return v if v else None
    return None


# ---------------------------------------------------------------------------
# Confidence normalisation
# ---------------------------------------------------------------------------

_CONFIDENCE_TYPE_BY_TOOL = {
    # FastTree's local-support values are SH-like, reported on [0, 1].
    "fasttree": "sh_like",
    # IQ-TREE with UFBoot (-B) produces ultrafast-bootstrap values on
    # [0, 100]. With -alrt only the values are SH-aLRT (also [0, 100]).
    # Repseq's default invokes UFBoot, so "ufboot" is the right label;
    # the writer is told the actual type via ``confidence_type=``.
    "iqtree": "ufboot",
}


def _confidence_type_for(tree_tool: str, override: Optional[str]) -> str:
    """Map a tree-builder name to its native support-value type.

    Accepts both lowercase keys (``"fasttree"``, ``"iqtree"``) and the
    display names the orchestrator passes for the ``<phylogeny>``
    element (``"FastTree"``, ``"IQ-TREE"``). Falls back to
    ``"bootstrap"`` if nothing matches.
    """
    if override and override != "auto":
        return override
    key = (tree_tool or "").lower().replace("-", "")
    return _CONFIDENCE_TYPE_BY_TOOL.get(key, "bootstrap")


def _normalize_confidence(raw: Optional[float], conf_type: str) -> Optional[int]:
    """Scale a tree-builder confidence value into a 0-100 integer.

    SH-like (FastTree) values come in [0, 1] and are rescaled.
    Everything else (UFBoot, SH-aLRT, classical bootstrap) is already
    on [0, 100]. A None input returns None — the caller skips the
    element rather than emitting a zero. Values are clamped to [0,
    100] in case a tool returns something slightly out of range.
    """
    if raw is None:
        return None
    if conf_type == "sh_like":
        # FastTree's [0,1] → [0,100]. Some FastTree builds emit
        # support on [0,100] already; we detect that by checking the
        # input range — anything >1 is left alone.
        value = float(raw) * 100.0 if 0.0 <= float(raw) <= 1.0 else float(raw)
    else:
        value = float(raw)
    value = max(0.0, min(100.0, value))
    return int(round(value))


# ---------------------------------------------------------------------------
# Element helpers
# ---------------------------------------------------------------------------

def _set_text(parent: ET.Element, tag: str, text: Optional[str]) -> Optional[ET.Element]:
    """Append a child element with text — but only if ``text`` is non-empty."""
    if text is None or str(text).strip() == "":
        return None
    elem = ET.SubElement(parent, tag)
    elem.text = str(text)
    return elem


def _format_branch_length(value: Optional[float]) -> Optional[str]:
    if value is None:
        return None
    # 6 significant figures is plenty for branch-length display and
    # avoids the long doubles Bio.Phylo sometimes produces.
    return f"{float(value):.6g}"


# ---------------------------------------------------------------------------
# Clade serialisation
# ---------------------------------------------------------------------------

def _collect_underlying_records(
    seq: Sequence,
) -> tuple[list[tuple[str, Optional[str]]], list[tuple[str, Optional[str], str]]]:
    """Return ``(nuc_records, protein_records)`` for one leaf sequence.

    * ``nuc_records``: ``[(accession, name)]`` per nucleotide segment,
      in segment order. For segmented CONCAT leaves the segments come
      from ``seq.concat_segments``; for non-segmented leaves the leaf
      itself is the single nuc record.
    * ``protein_records``: ``[(protein_id, product, source_segment_acc)]``
      flat across all segments, in (segment-order, CDS-order). The
      third tuple element lets the caller resolve which segment a
      protein lives on (currently unused but kept for future
      per-protein annotation).

    A leaf with no CDS info (``seq.proteins is None``) yields an empty
    protein_records list and the property emitter degrades gracefully.
    """
    segments = seq.concat_segments or [seq]
    nuc_records: list[tuple[str, Optional[str]]] = []
    protein_records: list[tuple[str, Optional[str], str]] = []
    for seg in segments:
        if seg.accession:
            nuc_records.append((seg.accession, _clean_name(seg.description or seg.header)))
        for p in (seg.proteins or []):
            pid = p.get("protein_id")
            if not pid:
                continue
            protein_records.append((pid, _clean_name(p.get("product")), seg.accession or ""))
    return nuc_records, protein_records


def _clean_name(value: Optional[str]) -> Optional[str]:
    """Trim a free-text name for emission inside ``<sequence><name>``.

    NCBI Virus FASTA headers store a pipe-separated metadata block
    right after the accession (``NC_078889.1 |species|host|...|segment``),
    so the parsed ``description`` field often starts with the
    separator. Stripping that here keeps the displayed name from
    starting with a stray ``|``. We also strip trailing pipes and
    collapse adjacent separators that result from empty metadata
    fields (``"...||..."`` → ``"... | ..."``).
    """
    if value is None:
        return None
    s = str(value).strip().strip("|").strip()
    # Collapse multi-pipe runs that came from empty NCBI-Virus fields.
    while "||" in s:
        s = s.replace("||", "|")
    return s or None


def _emit_sequence_element(
    clade: ET.Element,
    seq_type: str,
    accession: Optional[str],
    name: Optional[str],
) -> None:
    """Append one ``<sequence type="...">`` block to ``clade``.

    Skipped entirely when there's no accession AND no name, so we
    never write a bare ``<sequence/>`` that some viewers complain
    about.
    """
    if not accession and not name:
        return
    sequence_el = ET.SubElement(clade, "sequence", {"type": seq_type})
    if accession:
        acc_el = ET.SubElement(sequence_el, "accession", {"source": "ncbi"})
        acc_el.text = str(accession)
    _set_text(sequence_el, "name", name)


def _serialise_leaf(
    parent: ET.Element,
    leaf,
    seq: Sequence,
    label: str,
    conf_type: str,
) -> None:
    """Emit one terminal ``<clade>`` with rich annotation.

    The clade gets:

    * One ``<sequence>`` per underlying nucleotide segment AND one per
      protein CDS. Marker proteins (the ones used to build the tree)
      come first, then non-marker proteins, then the nucleotide
      segments. This ordering matters because many phyloXML renderers
      display only the first ``<sequence>`` element — putting the
      marker first means the visible label is the protein that
      actually drove the tree.
    * Three ``repseq:`` properties carrying comma-separated lists
      that reproduce the same data in a form renderers parse
      uniformly: ``repseq:nuc_acc``, ``repseq:protein_acc``, and
      ``repseq:protein_names``.
    * The isolate-level ``repseq:`` properties (host, country, date,
      strain, …) as before.
    """
    clade = ET.SubElement(parent, "clade")

    # PhyloXML schema order:
    # name -> branch_length -> confidence -> taxonomy -> sequence -> property
    _set_text(clade, "name", label)
    _set_text(clade, "branch_length", _format_branch_length(leaf.branch_length))

    conf_value = _normalize_confidence(getattr(leaf, "confidence", None), conf_type)
    if conf_value is not None:
        conf_el = ET.SubElement(clade, "confidence", {"type": conf_type})
        conf_el.text = str(conf_value)

    # <taxonomy>
    species = (
        seq.taxonomy.get_rank("species") if seq.taxonomy else None
    ) or seq.organism
    taxid = seq.taxonomy.taxid if seq.taxonomy else None
    if species or taxid is not None:
        tax_el = ET.SubElement(clade, "taxonomy")
        if taxid is not None:
            id_el = ET.SubElement(tax_el, "id", {"provider": "ncbi"})
            id_el.text = str(taxid)
        _set_text(tax_el, "scientific_name", species)

    # <sequence> elements: markers first, then other proteins, then nucs.
    nuc_records, protein_records = _collect_underlying_records(seq)
    marker_ids = list(seq.marker_protein_ids or [])
    marker_set = set(marker_ids)
    # Stable sort keeps marker ids in the configured marker order
    # (matters for segmented runs — L's marker before M's before S's).
    markers = [
        rec for pid in marker_ids
        for rec in protein_records if rec[0] == pid
    ]
    others = [rec for rec in protein_records if rec[0] not in marker_set]
    for pid, product, _src in markers:
        _emit_sequence_element(clade, "protein", pid, product)
    for pid, product, _src in others:
        _emit_sequence_element(clade, "protein", pid, product)
    for acc, name in nuc_records:
        _emit_sequence_element(clade, "dna", acc, name)

    # repseq: summary properties — same data as the <sequence>
    # elements, comma-joined, for renderers that only show the first
    # <sequence>. nuc_acc lists segments in segment order; protein_acc
    # / protein_names follow the same ordering as the <sequence>
    # blocks above (markers first, then others).
    ordered_protein_ids = [pid for pid, _p, _s in markers + others]
    ordered_protein_names = [
        product for _pid, product, _s in markers + others if product
    ]
    ordered_nuc_accs = [acc for acc, _n in nuc_records]
    for ref, values in (
        ("repseq:nuc_acc", ordered_nuc_accs),
        ("repseq:protein_acc", ordered_protein_ids),
        ("repseq:protein_names", ordered_protein_names),
    ):
        if not values:
            continue
        prop = ET.SubElement(
            clade,
            "property",
            {"ref": ref, "datatype": "xsd:string", "applies_to": "clade"},
        )
        prop.text = ", ".join(values)

    # Isolate-level <property> elements (host, country, date, …).
    for ref, key in _LEAF_PROPERTIES:
        value = _leaf_property_value(seq, key)
        if not value:
            continue
        prop = ET.SubElement(
            clade,
            "property",
            {
                "ref": ref,
                "datatype": "xsd:string",
                "applies_to": "clade",
            },
        )
        prop.text = str(value)


def _serialise_internal(
    parent: ET.Element,
    node,
    conf_type: str,
    seq_by_id: dict[str, Sequence],
    label_by_id: dict[str, str],
) -> None:
    """Emit one internal ``<clade>`` and recurse into children."""
    clade = ET.SubElement(parent, "clade")

    # PhyloXML schema order on internals:
    #   name -> branch_length -> confidence -> taxonomy -> child clades
    # Both <name> and <taxonomy> are only emitted when the LCA
    # annotator (repseq.phylo.lca) has set _lca_name on the node.
    lca_name = getattr(node, "_lca_name", None)
    if lca_name:
        _set_text(clade, "name", lca_name)
    _set_text(clade, "branch_length", _format_branch_length(node.branch_length))
    conf_value = _normalize_confidence(getattr(node, "confidence", None), conf_type)
    if conf_value is not None:
        conf_el = ET.SubElement(clade, "confidence", {"type": conf_type})
        conf_el.text = str(conf_value)
    if lca_name:
        tax_el = ET.SubElement(clade, "taxonomy")
        _set_text(tax_el, "scientific_name", lca_name)
        _set_text(tax_el, "rank", phyloxml_rank(getattr(node, "_lca_rank", None)))

    for child in node.clades:
        if child.is_terminal():
            original = (
                getattr(child, "_repseq_original_id", None)
                or child.name
            )
            seq = seq_by_id.get(original)
            if seq is None:
                # Fall back to a minimal leaf — should not happen if
                # id_map is correct, but better safe than crash.
                fallback = ET.SubElement(clade, "clade")
                _set_text(fallback, "name", original or "")
                _set_text(
                    fallback,
                    "branch_length",
                    _format_branch_length(child.branch_length),
                )
                continue
            _serialise_leaf(
                clade,
                child,
                seq,
                label=label_by_id.get(original, original),
                conf_type=conf_type,
            )
        else:
            _serialise_internal(
                clade,
                child,
                conf_type=conf_type,
                seq_by_id=seq_by_id,
                label_by_id=label_by_id,
            )


# ---------------------------------------------------------------------------
# Phylogeny-level header
# ---------------------------------------------------------------------------

def _build_phylogeny_name(
    prefix: str,
    alphabet: str,
    msa_tool: str,
    tree_tool: str,
    model: Optional[str],
) -> str:
    """Format the ``<phylogeny><name>`` element.

    Shape: ``"{prefix} [{alphabet}|{msa_tool}|{tree_tool} {model}]"``.
    """
    parts = [tree_tool]
    if model:
        parts.append(str(model))
    return f"{prefix} [{alphabet}|{msa_tool}|{' '.join(parts)}]"


def _markers_summary(reps: list[Sequence]) -> Optional[str]:
    """Summarise the marker proteins that fed the tree, for the
    ``<phylogeny><description>`` element.

    Two shapes, depending on input:

    * **Segmented**: ``"L:polymerase, M:glycoprotein, S:nucleoprotein"``.
      Segment order is taken from the first rep's ``concat_segments``
      (the same order across every isolate, set by
      ``build_concatenated_sequences`` from the virus config). When
      different isolates picked different products on the same
      segment, those products are pipe-joined and alphabetised within
      the segment: ``"L:polymerase|RdRp"``.
    * **Non-segmented**: ``"polymerase, NS5"`` — the unique product
      names selected as marker across all reps, alphabetised.

    Returns ``None`` when no rep has ``marker_protein_ids`` set
    (nucleotide-alphabet runs, or protein runs where every rep
    fell through to no-marker). The caller then omits the
    ``markers=...`` field from the description.
    """
    has_concat = any(r.concat_segments for r in reps)

    if has_concat:
        per_segment: dict[str, set[str]] = {}
        segment_order: list[str] = []
        seen_segments: set[str] = set()
        for rep in reps:
            if not rep.concat_segments or not rep.marker_protein_ids:
                continue
            marker_set = set(rep.marker_protein_ids)
            for seg in rep.concat_segments:
                if not seg.segment:
                    continue
                if seg.segment not in seen_segments:
                    segment_order.append(seg.segment)
                    seen_segments.add(seg.segment)
                for protein in (seg.proteins or []):
                    if protein.get("protein_id") in marker_set:
                        product = protein.get("product")
                        if product:
                            per_segment.setdefault(seg.segment, set()).add(product)
        bits: list[str] = []
        for seg_name in segment_order:
            products = per_segment.get(seg_name)
            if not products:
                continue
            bits.append(f"{seg_name}:{'|'.join(sorted(products))}")
        return ", ".join(bits) if bits else None

    products: set[str] = set()
    for rep in reps:
        if not rep.marker_protein_ids or not rep.proteins:
            continue
        marker_set = set(rep.marker_protein_ids)
        for protein in rep.proteins:
            if protein.get("protein_id") in marker_set:
                p = protein.get("product")
                if p:
                    products.add(p)
    return ", ".join(sorted(products)) if products else None


def _build_phylogeny_description(
    alphabet: str,
    msa_tool: str,
    msa_version: str,
    tree_tool: str,
    tree_version: str,
    model: Optional[str],
    ufboot: Optional[int],
    extra_msa_args: list[str],
    extra_tree_args: list[str],
    rooting_method: Optional[str] = None,
    markers: Optional[str] = None,
) -> str:
    """Compose the ``<phylogeny><description>`` element."""
    when = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    bits = [
        f"repseq {REPSEQ_VERSION} phylogeny",
        f"generated {when}",
        f"alphabet={alphabet}",
        f"MSA={msa_tool} {msa_version}",
    ]
    if extra_msa_args:
        bits.append(f"mafft_args={' '.join(extra_msa_args)}")
    bits.append(f"tree={tree_tool} {tree_version}")
    if model:
        bits.append(f"model={model}")
    bits.append(f"bootstrap={ufboot if ufboot else 'none'}")
    if extra_tree_args:
        bits.append(f"tree_args={' '.join(extra_tree_args)}")
    if rooting_method:
        bits.append(f"rooting={rooting_method}")
    if markers:
        bits.append(f"markers={markers}")
    return " | ".join(bits)


# ---------------------------------------------------------------------------
# Top-level write
# ---------------------------------------------------------------------------

def write_phyloxml(
    newick_path: Optional[Path],
    phyloxml_path: Path,
    representatives: list[Sequence],
    id_map: dict[str, str],
    *,
    cfg: dict[str, Any],
    prefix: str,
    alphabet: str,
    msa_tool: str,
    msa_version: str,
    tree_tool: str,
    tree_version: str,
    model: Optional[str],
    ufboot: Optional[int],
    extra_msa_args: Optional[list[str]] = None,
    extra_tree_args: Optional[list[str]] = None,
    tree=None,
    rooting_method: Optional[str] = None,
) -> None:
    """Render a tree to a richly-annotated phyloXML file at
    ``phyloxml_path``.

    Either ``newick_path`` or a pre-parsed ``tree`` must be provided —
    the caller (the pipeline orchestrator) parses + roots + LCA-annotates
    the tree first, then hands it in here for serialisation. Falling
    back to parsing a Newick path directly is supported for callers
    that don't need rooting/LCA.

    ``representatives`` is the list of sequences whose ids appear as
    leaves; ``id_map`` maps the short ids (S0001…) back to the
    corresponding ``seq.id`` values. The writer ladderizes the tree
    (reverse=True, larger clades on top) before emitting XML.
    ``rooting_method`` (when provided) is included in the
    ``<phylogeny><description>`` so the user can see which method
    actually fired in the auto chain.
    """
    if extra_msa_args is None:
        extra_msa_args = []
    if extra_tree_args is None:
        extra_tree_args = []

    if tree is None:
        if newick_path is None:
            raise ValueError("write_phyloxml needs either tree or newick_path")
        tree = Phylo.read(str(newick_path), "newick")
    tree.ladderize(reverse=True)

    # Map short_id (as found in the Newick) → original seq.id, and
    # record the original id on the terminal clade for downstream
    # serialisation. We don't rename the clade itself because some
    # original ids contain pipes — keeping the original on a private
    # attribute lets the leaf serialiser format its own label.
    for terminal in tree.get_terminals():
        if terminal.name in id_map:
            terminal._repseq_original_id = id_map[terminal.name]

    segmented = bool((cfg or {}).get("segmented", {}).get("enabled"))
    label_format = pick_format_string(cfg, segmented=segmented)
    label_opts = labeling_options(cfg)
    label_by_id = {
        seq.id: format_leaf_label(seq, label_format, **label_opts)
        for seq in representatives
    }
    seq_by_id = {seq.id: seq for seq in representatives}

    phyloxml_cfg = (cfg or {}).get("phylo", {}).get("phyloxml", {}) or {}
    conf_type = _confidence_type_for(
        tree_tool, phyloxml_cfg.get("confidence_type"),
    )

    root_attrs = {
        "xmlns": PHYLOXML_NS,
        "xmlns:xsi": XSI_NS,
        "xmlns:repseq": REPSEQ_NS,
        "xsi:schemaLocation": (
            f"{PHYLOXML_NS} http://www.phyloxml.org/1.10/phyloxml.xsd"
        ),
    }
    root = ET.Element("phyloxml", root_attrs)

    phylogeny = ET.SubElement(
        root,
        "phylogeny",
        {"rooted": "true" if getattr(tree, "rooted", False) else "false"},
    )
    _set_text(
        phylogeny,
        "name",
        _build_phylogeny_name(prefix, alphabet, msa_tool, tree_tool, model),
    )
    _set_text(
        phylogeny,
        "description",
        _build_phylogeny_description(
            alphabet=alphabet,
            msa_tool=msa_tool,
            msa_version=msa_version,
            tree_tool=tree_tool,
            tree_version=tree_version,
            model=model,
            ufboot=ufboot,
            extra_msa_args=extra_msa_args,
            extra_tree_args=extra_tree_args,
            rooting_method=rooting_method,
            markers=_markers_summary(representatives),
        ),
    )

    # Walk from the tree root. Bio.Phylo's root is itself a clade; we
    # treat it as an "internal" clade and let serialisation recurse.
    _serialise_internal(
        phylogeny,
        tree.root,
        conf_type=conf_type,
        seq_by_id=seq_by_id,
        label_by_id=label_by_id,
    )

    # Pretty-print: stdlib doesn't ship a pretty printer that handles
    # namespaces correctly, so we run a small one. Indents are 2
    # spaces, mirroring most phyloXML examples in the wild.
    _indent(root)

    phyloxml_path.parent.mkdir(parents=True, exist_ok=True)
    tree_doc = ET.ElementTree(root)
    tree_doc.write(
        str(phyloxml_path),
        encoding="utf-8",
        xml_declaration=True,
    )


def _indent(elem: ET.Element, level: int = 0) -> None:
    """In-place pretty-print for an ElementTree (stdlib has no equivalent
    that round-trips namespaces, so we re-implement)."""
    pad = "\n" + "  " * level
    child_pad = "\n" + "  " * (level + 1)
    if len(elem):
        if not (elem.text or "").strip():
            elem.text = child_pad
        for i, child in enumerate(elem):
            _indent(child, level + 1)
            if not (child.tail or "").strip():
                child.tail = child_pad if i < len(elem) - 1 else pad
    else:
        if level and not (elem.tail or "").strip():
            elem.tail = pad
