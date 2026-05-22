"""User-facing exceptions that map to friendly, traceback-free CLI messages.

Anything raised as a :class:`RepseqError` (or subclass) is caught at the CLI
boundary (``repseq.cli._RepseqGroup.main``), printed as a single
``Error: ...`` line to stderr, and exits non-zero **without** a Python
traceback. These are reserved for mistakes a bench scientist can actually
fix — a wrong path, a malformed config file, a missing external tool.

Genuine internal bugs are deliberately *not* funnelled through here: they
propagate as a normal traceback so they get noticed and reported. The
boundary additionally renders a small set of known external-tool failures
(MMseqs2 / cd-hit not on PATH) as friendly errors even though those tool
classes live outside this hierarchy — see ``cli._RepseqGroup``.
"""

from __future__ import annotations


class RepseqError(Exception):
    """Base for user-facing errors rendered without a traceback.

    The string form is shown verbatim after an ``Error: `` prefix, so the
    message should be plain English and name the next step. Multi-line
    messages are fine; indent continuation lines with 7 spaces so they
    align under the prefix.
    """


class ConfigError(RepseqError):
    """A problem with the config file: missing, unparseable, or wrong shape."""


class InputError(RepseqError):
    """A problem with an input FASTA file: missing, unreadable, empty, or
    not actually FASTA."""
