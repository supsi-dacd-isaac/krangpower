import os
import sys
from contextlib import contextmanager
import platform


@contextmanager
def stdout_redirected(to=os.devnull):
    """
    import os

    with stdout_redirected(to=filename):
        print("from Python")
        os.system("echo non-Python applications are also supported")
    """
    if platform.system() == 'Windows':
        pass  # nothing like this works on win

    else:
        fd = sys.stdout.fileno()

        def _redirect_stdout(to):
            sys.stdout.close()  # + implicit flush()
            os.dup2(to.fileno(), fd)  # fd writes to 'to' file
            sys.stdout = os.fdopen(fd, 'w')  # Python writes to fd

        with os.fdopen(os.dup(fd), 'w') as old_stdout:
            with open(to, 'w') as file:
                _redirect_stdout(to=file)

            try:
                yield  # allow code to be run with the redirected stdout
            finally:
                _redirect_stdout(to=old_stdout)  # restore stdout. buffering and flags such as CLOEXEC may be different
