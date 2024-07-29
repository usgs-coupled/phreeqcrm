import pytest
import os
import shutil

from constants import FilePaths

def pytest_sessionstart():
    cwd = os.getcwd()

    filenames = [FilePaths.DATABASE, FilePaths.PQI]
    for filename in filenames:
        src = filename
        dest = os.path.join(cwd, os.path.basename(filename))
        if not os.path.exists(dest):
            shutil.copy(src, dest)
