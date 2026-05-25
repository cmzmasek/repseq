"""run_streaming — tee stderr live AND buffer it for error reporting.

The real subprocesses (MAFFT/IQ-TREE/FastTree) are mocked at the
``run_mafft``/``run_iqtree``/``run_fasttree`` level in the rest of the
suite, so this file is the only place we exercise the streaming helper
directly. Tests use ``python -c "..."`` as the child to keep them
portable, hermetic, and fast.
"""
from __future__ import annotations

import io
import sys
from pathlib import Path

import pytest

from repseq.utils.subprocess_stream import (
    StreamedProcessError,
    run_streaming,
)


def test_run_streaming_returns_buffered_stderr():
    cmd = [sys.executable, "-c", "import sys; sys.stderr.write('hello\\n')"]
    sink = io.StringIO()
    out = run_streaming(cmd, stderr_dest=sink)
    assert "hello" in out
    assert "hello" in sink.getvalue()


def test_run_streaming_writes_stream_prefix_per_line():
    cmd = [
        sys.executable, "-c",
        "import sys; sys.stderr.write('a\\nb\\n')",
    ]
    sink = io.StringIO()
    run_streaming(cmd, stream_prefix="[tag] ", stderr_dest=sink)
    streamed = sink.getvalue()
    # Each source line picks up the prefix.
    assert "[tag] a\n" in streamed
    assert "[tag] b\n" in streamed


def test_run_streaming_routes_stdout_to_file(tmp_path):
    out_file = tmp_path / "out.txt"
    cmd = [sys.executable, "-c", "print('PAYLOAD')"]
    run_streaming(cmd, stdout_file=out_file)
    assert out_file.read_text().strip() == "PAYLOAD"


def test_run_streaming_raises_on_nonzero_exit_with_buffered_stderr():
    cmd = [
        sys.executable, "-c",
        "import sys; sys.stderr.write('boom\\n'); sys.exit(7)",
    ]
    sink = io.StringIO()
    with pytest.raises(StreamedProcessError) as exc:
        run_streaming(cmd, stderr_dest=sink)
    assert exc.value.returncode == 7
    assert "boom" in (exc.value.stderr or "")
    # And the stream sink got the same content (the user saw it live).
    assert "boom" in sink.getvalue()


def test_run_streaming_check_false_returns_instead_of_raising():
    cmd = [
        sys.executable, "-c",
        "import sys; sys.stderr.write('warn\\n'); sys.exit(2)",
    ]
    sink = io.StringIO()
    out = run_streaming(cmd, check=False, stderr_dest=sink)
    assert "warn" in out


def test_run_streaming_creates_stdout_file_parent_dir(tmp_path):
    """A leading non-existent dir component shouldn't blow up — the
    helper mkdir(parents=True)'s the parent before opening."""
    target = tmp_path / "nested" / "deeper" / "out.txt"
    cmd = [sys.executable, "-c", "print('ok')"]
    run_streaming(cmd, stdout_file=target)
    assert target.read_text().strip() == "ok"
