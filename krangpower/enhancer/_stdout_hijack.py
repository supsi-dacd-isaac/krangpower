import os
import sys
from platform import system
from tempfile import TemporaryFile
from threading import Thread

from .._config_loader import ODSS_STDOUT_SUPPRESSED


class NullCm:
    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


class _stdout_redirected_linux:

    # the point of this class is the temporary redirection of stdout to a pipe that's continuously drained
    # by a thread. It's clean, straightforward and it doesn't work under Windows.

    @staticmethod
    def drain_pipe(pipe_read_side):
        while True:
            data = os.read(pipe_read_side, 1024)
            if not data:
                break

    def __enter__(self):
        self.stdout_fileno = sys.stdout.fileno()
        self.stdout_save = os.dup(self.stdout_fileno)  # refers to the pre-existent stdout file
        self.stdout_pipe = os.pipe()
        os.dup2(self.stdout_pipe[1], self.stdout_fileno)  # hijacks stdout fileno's handler to the pipe's write end
        os.close(self.stdout_pipe[1])  # we need no more the open pipe's write end

        self.t = Thread(target=self.drain_pipe, args=(self.stdout_pipe[0],))
        self.t.start()

    def __exit__(self, exc_type, exc_val, exc_tb):
        os.close(self.stdout_fileno)  # this triggers the thread's exit condition
        self.t.join()  # waits for the thread to exit

        # clean up
        os.close(self.stdout_pipe[0])  # pipe's read end
        os.dup2(self.stdout_save, self.stdout_fileno)  # resets stdout to the original file handler
        os.close(self.stdout_save)


class _stdout_redirected_win:

    def __enter__(self):
        self._org = sys.stdout
        sys.stdout = sys.__stdout__
        fdout = sys.stdout.fileno()
        self._file = TemporaryFile()
        self._dup = None
        if fdout >= 0:
            self._dup = os.dup(fdout)
            os.dup2(self._file.fileno(), fdout)

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.flush()
        if self._dup is not None:
            os.dup2(self._dup, sys.stdout.fileno())
            os.close(self._dup)
        sys.stdout = self._org
        self._file.seek(0)
        self._file.close()


# when loading the module, the exposed context manager is chosen accordingly
if not ODSS_STDOUT_SUPPRESSED:
    stdout_redirected = NullCm

elif system() == 'Windows':
    stdout_redirected = _stdout_redirected_win

elif system() == 'Linux':
    stdout_redirected = _stdout_redirected_linux

else:
    raise OSError('Unsupported system: {}'.format(system()))
