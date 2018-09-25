import os
import sys
import tempfile
import threading
from platform import system


class stdout_redirected_linux:

    @staticmethod
    def drain_pipe(pipe_read_side):
        while True:
            data = os.read(pipe_read_side, 1024)
            if not data:
                break

    def __enter__(self):
        self.stdout_fileno = sys.stdout.fileno()
        self.stdout_save = os.dup(stdout_fileno)
        self.stdout_pipe = os.pipe()
        os.dup2(self.stdout_pipe[1], self.stdout_fileno)
        os.close(self.stdout_pipe[1])

        self.t = threading.Thread(target=drain_pipe, args=(self.stdout_pipe[0],))
        self.t.start()

    def __exit__(self, exc_type, exc_val, exc_tb):
        os.close(self.stdout_fileno)
        self.t.join()
        os.close(self.stdout_pipe[0])
        os.dup2(self.stdout_save, self.stdout_fileno)
        os.close(self.stdout_save)


class stdout_redirected_win:

    def __enter__(self):
        self._org = sys.stdout
        sys.stdout = sys.__stdout__
        fdout = sys.stdout.fileno()
        self._file = tempfile.TemporaryFile()
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


if system() == 'Windows':
    stdout_redirected = stdout_redirected_win

elif system() == 'Linux':
    stdout_redirected = stdout_redirected_linux

else:
    raise OSError('Unsupported system: {}'.format(system()))
