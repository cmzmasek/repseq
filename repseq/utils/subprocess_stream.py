"""Run a subprocess while streaming its stderr live to the user.

The pre-existing pattern in this repo was
``subprocess.run(..., capture_output=True, check=True)`` — clean and
mockable, but on a long-running step (MAFFT on thousands of sequences,
IQ-TREE ModelFinder, a slow hmmscan) the terminal goes silent for
minutes at a time, and a bench scientist reasonably wonders whether the
pipeline froze. This helper splits the difference: stderr lines stream
to the user's terminal as they arrive, AND each line is buffered so the
caller can include the full text in an error message if the subprocess
exits non-zero (the behaviour callers relied on).

Used by the MAFFT / IQ-TREE / FastTree wrappers. Tests mock these
wrappers at a higher level, so the streaming machinery itself is
covered by a small set of direct tests in ``tests/test_subprocess_stream.py``.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Optional


class StreamedProcessError(subprocess.CalledProcessError):
    """Same shape as :class:`subprocess.CalledProcessError`. Exists so the
    helper's failure type is unambiguous in tracebacks even though we
    inherit all its attributes from the parent class."""


def run_streaming(
    argv: list[str],
    *,
    stdout_file: Optional[Path] = None,
    cwd: Optional[Path] = None,
    check: bool = True,
    stream_prefix: str = "",
    stderr_dest=None,
    stream_stderr: bool = True,
) -> str:
    """Run ``argv`` and stream its stderr live; buffer for error reporting.

    Args:
        argv: command + arguments.
        stdout_file: when given, the child's stdout is written there (used
            by MAFFT/FastTree which print their primary output to stdout).
            When ``None``, stdout is discarded.
        cwd: working directory for the child.
        check: when True (default), raise :class:`StreamedProcessError`
            on a non-zero exit. When False, return normally and the
            caller checks the buffered stderr.
        stream_prefix: prepended to every streamed line (e.g.
            ``"[mafft] "``) so a glance at the log identifies the source
            when multiple subprocesses run in sequence.
        stderr_dest: file-like object the live stream is written to.
            Defaults to ``sys.stderr``; tests pass a ``StringIO`` to
            assert on the streamed output.
        stream_stderr: when True (default), each stderr line is echoed
            to ``stderr_dest`` as it arrives (the heartbeat). When
            False, stderr is still consumed line-by-line and buffered
            so the on-failure error message keeps its full text, but
            nothing is written to ``stderr_dest`` during the run — the
            terminal stays quiet. The ``--verbose`` CLI flag toggles
            this in the MAFFT / IQ-TREE / FastTree wrappers.

    Returns:
        The full buffered stderr text as a single string.

    Raises:
        StreamedProcessError: when ``check=True`` and the child exits
            non-zero. The exception's ``stderr`` attribute is the same
            buffered text the function would return.
    """
    if stderr_dest is None:
        stderr_dest = sys.stderr

    if stdout_file is not None:
        stdout_file = Path(stdout_file)
        stdout_file.parent.mkdir(parents=True, exist_ok=True)
        stdout_handle = open(stdout_file, "w")
        stdout_target = stdout_handle
    else:
        stdout_handle = None
        stdout_target = subprocess.DEVNULL

    buf: list[str] = []
    try:
        proc = subprocess.Popen(
            argv,
            stdout=stdout_target,
            stderr=subprocess.PIPE,
            cwd=str(cwd) if cwd else None,
            text=True,
            bufsize=1,  # line-buffered
        )
        assert proc.stderr is not None  # subprocess.PIPE guarantees this
        for line in proc.stderr:
            buf.append(line)
            if stream_stderr:
                if stream_prefix:
                    stderr_dest.write(stream_prefix + line)
                else:
                    stderr_dest.write(line)
                try:
                    stderr_dest.flush()
                except Exception:
                    # A test StringIO or a closed terminal must not break the run.
                    pass
        rc = proc.wait()
    finally:
        if stdout_handle is not None:
            stdout_handle.close()

    full_stderr = "".join(buf)
    if check and rc != 0:
        raise StreamedProcessError(rc, argv, output=None, stderr=full_stderr)
    return full_stderr
