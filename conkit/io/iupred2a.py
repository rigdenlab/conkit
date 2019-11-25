import pandas as pd

class DisorderFile(object):

    """Class to contain disorder prediction data from a file"""

    def __init__(self):
        self._fname = None
        self._iupred2a = None

    @property
    def iupred2a(self):
        return self._iupred2a

    @iupred2a.setter
    def iupred2a(self, value):
        self._iupred2a = value

    @property
    def fname(self):
        return self._fname

    @fname.setter
    def fname(self, value):
        self._fname = value


class iupred2aParser(object):
    """Class for parsing iupred2a files."""

    def __init__(self):
        """Initialise an instance of the iupred2aParser."""

        pass

    @staticmethod
    def read(fname):
        """Return iupred2a prediction instance."""
        iupred_data = pd.read_csv(fname, sep='\t', comment='#', skip_blank_lines=True, header=None)
        iupred2a = iupred_data[2].tolist()

        hierarchy = DisorderFile()
        hierarchy.fname = fname
        hierarchy.iupred2a = iupred2a
        return hierarchy

    @staticmethod
    def write(self):
        raise NotImplementedError('This is not implemented!')


if __name__ == '__main__':

    topcons_prediction = iupred2aParser().read('/Users/shahrammesdaghi/cmap_plus/iupred2a.txt')
    print(topcons_prediction.iupred2a)
