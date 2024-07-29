import pytest
import os

from constants import FilePaths

def pytest_sessionstart():
    cwd = os.getcwd()
    print(f"pytest_sessionstart -- Current working directory: {cwd}")

    filenames = [FilePaths.DATABASE, FilePaths.PQI]

    for filename in filenames:
        src = filename
        dest = os.path.join(cwd, os.path.basename(filename))
        if not os.path.exists(dest):
            shutil.copy(src, dest)
        else:
            print(f"pytest_sessionstart -- Destination already exists: {dest}")
