"""Persistent SQLite cache for taxonomy and metadata lookups."""

from __future__ import annotations

import json
import sqlite3
import threading
import time
from pathlib import Path
from typing import Any, Optional


class TaxonomyCache:
    """Thread-safe SQLite-backed cache shared across runs.

    Keys are (source, accession_or_taxid) pairs.
    Values are JSON-serialised dicts.
    """

    def __init__(self, cache_dir: str | Path, ttl_days: int = 90) -> None:
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(parents=True, exist_ok=True)
        self._db_path = cache_dir / "taxonomy.db"
        self._ttl_seconds = ttl_days * 86400
        self._conn = sqlite3.connect(str(self._db_path), check_same_thread=False)
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._lock = threading.Lock()
        self._create_tables()

    def _create_tables(self) -> None:
        self._conn.execute(
            """
            CREATE TABLE IF NOT EXISTS cache (
                source      TEXT NOT NULL,
                key         TEXT NOT NULL,
                value       TEXT NOT NULL,
                cached_at   INTEGER NOT NULL,
                PRIMARY KEY (source, key)
            )
            """
        )
        self._conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_cache_source ON cache(source)"
        )
        self._conn.commit()

    # ------------------------------------------------------------------
    # CRUD
    # ------------------------------------------------------------------

    def get(self, source: str, key: str) -> Optional[dict[str, Any]]:
        """Return cached value or None if missing/expired."""
        with self._lock:
            row = self._conn.execute(
                "SELECT value, cached_at FROM cache WHERE source=? AND key=?",
                (source, key),
            ).fetchone()
            if row is None:
                return None
            value_json, cached_at = row
            if self._ttl_seconds > 0 and (time.time() - cached_at) > self._ttl_seconds:
                self._conn.execute(
                    "DELETE FROM cache WHERE source=? AND key=?", (source, key)
                )
                self._conn.commit()
                return None
        return json.loads(value_json)

    def set(self, source: str, key: str, value: dict[str, Any]) -> None:
        with self._lock:
            self._conn.execute(
                """
                INSERT OR REPLACE INTO cache (source, key, value, cached_at)
                VALUES (?, ?, ?, ?)
                """,
                (source, key, json.dumps(value), int(time.time())),
            )
            self._conn.commit()

    def delete(self, source: str, key: str) -> None:
        with self._lock:
            self._conn.execute(
                "DELETE FROM cache WHERE source=? AND key=?", (source, key)
            )
            self._conn.commit()

    # ------------------------------------------------------------------
    # Management
    # ------------------------------------------------------------------

    def clear(self, source: Optional[str] = None) -> int:
        """Delete all entries, optionally filtered by source. Returns rows deleted."""
        with self._lock:
            if source:
                cur = self._conn.execute(
                    "DELETE FROM cache WHERE source=?", (source,)
                )
            else:
                cur = self._conn.execute("DELETE FROM cache")
            self._conn.commit()
            return cur.rowcount

    def stats(self) -> dict[str, Any]:
        """Return cache statistics."""
        with self._lock:
            total = self._conn.execute("SELECT COUNT(*) FROM cache").fetchone()[0]
            by_source = {
                row[0]: row[1]
                for row in self._conn.execute(
                    "SELECT source, COUNT(*) FROM cache GROUP BY source"
                ).fetchall()
            }
        size_bytes = self._db_path.stat().st_size if self._db_path.exists() else 0
        return {
            "total_entries": total,
            "by_source": by_source,
            "db_size_mb": round(size_bytes / 1_048_576, 2),
            "db_path": str(self._db_path),
        }

    def purge_expired(self) -> int:
        """Remove expired entries. Returns rows deleted."""
        cutoff = int(time.time()) - self._ttl_seconds
        with self._lock:
            cur = self._conn.execute(
                "DELETE FROM cache WHERE cached_at < ?", (cutoff,)
            )
            self._conn.commit()
            return cur.rowcount

    def close(self) -> None:
        self._conn.close()
