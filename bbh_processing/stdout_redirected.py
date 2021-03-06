"""
with stdout_redirected():
    CODE
Suppresses output from CODE.
"""

import os
import sys
import contextlib
def fileno(file_or_fd):
    """
    Return file number.
    """
    fds = getattr(file_or_fd, 'fileno', lambda: file_or_fd)()
    if not isinstance(fds, int):
        raise ValueError("Expected file or file descriptor")
    return fds

@contextlib.contextmanager
def stdout_redirected(too=os.devnull, stdout=None):
    """
    http://stackoverflow.com/a/22434262/190597 (J.F. Sebastian)
    """
    if stdout is None:
        stdout = sys.stdout

    stdout_fd = fileno(stdout)
    # copy stdout_fd before it is overwritten
    #NOTE: `copied` is inheritable on Windows when duplicating a standard stream
    with os.fdopen(os.dup(stdout_fd), 'wb') as copied:
        stdout.flush()  # flush library buffers that dup2 knows nothing about
        try:
            os.dup2(fileno(too), stdout_fd)  # $ exec >&to
        except ValueError:  # filename
            with open(too, 'wb') as to_file:
                os.dup2(to_file.fileno(), stdout_fd)  # $ exec > to
        try:
            yield stdout # allow code to be run with the redirected stdout
        finally:
            # restore stdout to its previous value
            #NOTE: dup2 makes stdout_fd inheritable unconditionally
            stdout.flush()
            os.dup2(copied.fileno(), stdout_fd)  # $ exec >&copiedi
