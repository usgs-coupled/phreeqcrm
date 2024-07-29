# constants.py

import os

class FilePaths:
    file_directory = os.path.dirname(__file__)
    YAML     = os.path.join(file_directory, 'AdvectBMI_py.yaml')
    DATABASE = os.path.join(file_directory, 'phreeqc.dat')
    PQI      = os.path.join(file_directory, 'advect.pqi')
