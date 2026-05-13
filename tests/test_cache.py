"""TaxonomyCache: round-trip, TTL expiry, clear/purge stats."""
from __future__ import annotations

import time

from repseq.taxonomy.cache import TaxonomyCache


def test_cache_set_get_roundtrip(tmp_cache_dir):
    cache = TaxonomyCache(tmp_cache_dir)
    cache.set("ncbi", "MW626064.1", {"organism": "Influenza A virus", "taxid": 11320})
    out = cache.get("ncbi", "MW626064.1")
    assert out == {"organism": "Influenza A virus", "taxid": 11320}


def test_cache_miss_returns_none(tmp_cache_dir):
    cache = TaxonomyCache(tmp_cache_dir)
    assert cache.get("ncbi", "nonexistent") is None


def test_cache_ttl_expiry(tmp_cache_dir, monkeypatch):
    """A TTL-expired entry should be returned as None and removed."""
    cache = TaxonomyCache(tmp_cache_dir, ttl_days=1)

    # Stamp the entry as if cached two days ago.
    cache._conn.execute(
        "INSERT OR REPLACE INTO cache (source, key, value, cached_at) VALUES (?, ?, ?, ?)",
        ("ncbi", "k1", '{"a": 1}', int(time.time()) - 2 * 86400),
    )
    cache._conn.commit()

    assert cache.get("ncbi", "k1") is None
    # And it should have been removed
    row = cache._conn.execute(
        "SELECT value FROM cache WHERE source=? AND key=?", ("ncbi", "k1")
    ).fetchone()
    assert row is None


def test_cache_clear_by_source(tmp_cache_dir):
    cache = TaxonomyCache(tmp_cache_dir)
    cache.set("ncbi", "a", {"v": 1})
    cache.set("uniprot", "b", {"v": 2})
    deleted = cache.clear("ncbi")
    assert deleted == 1
    assert cache.get("ncbi", "a") is None
    assert cache.get("uniprot", "b") == {"v": 2}


def test_cache_purge_expired(tmp_cache_dir):
    cache = TaxonomyCache(tmp_cache_dir, ttl_days=1)
    # one fresh entry, one expired
    cache.set("ncbi", "fresh", {"v": 1})
    cache._conn.execute(
        "INSERT OR REPLACE INTO cache (source, key, value, cached_at) VALUES (?, ?, ?, ?)",
        ("ncbi", "stale", '{"v": 2}', int(time.time()) - 2 * 86400),
    )
    cache._conn.commit()

    n = cache.purge_expired()
    assert n == 1
    assert cache.get("ncbi", "fresh") == {"v": 1}


def test_cache_stats_counts_per_source(tmp_cache_dir):
    cache = TaxonomyCache(tmp_cache_dir)
    cache.set("ncbi", "a", {"v": 1})
    cache.set("ncbi", "b", {"v": 2})
    cache.set("uniprot", "c", {"v": 3})
    stats = cache.stats()
    assert stats["total_entries"] == 3
    assert stats["by_source"] == {"ncbi": 2, "uniprot": 1}
