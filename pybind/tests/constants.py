# constants.py
class Constants:
    def __init__(self):
        file_directory = os.path.dirname(__file__)
        self._yaml     = os.path.join(file_directory, 'AdvectBMI_py.yaml')
        self._database = os.path.join(file_directory, 'phreeqc.dat')
        self._pqi      = os.path.join(file_directory, 'advect_pqi')

    @property
    def yaml(self):
        return self._yaml

    @property
    def database(self):
        return self._database

    @property
    def pqi(self):
        return self._pqi
