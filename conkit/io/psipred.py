import pandas as pd


class SspredFile(object):

    """Class to contain SS prediction data from a file"""

    def __init__(self):
        self._fname = None
        self._df = None

    @property
    def df(self):
        return self._df

    @df.setter
    def df(self, value):
        self._df = value

    @property
    def fname(self):
        return self._fname

    @fname.setter
    def fname(self, value):
        self._fname = value

    def __getitem__(self, item):
        return self.df[item].tolist()


class PsipredParser(object):

    """Class for parsing psipred files."""

    def __init__(self):

        """Initialise an instance of the psipredParser."""

        pass

    @staticmethod
    def read(fname):
        """Return psipred prediction instance."""

        ss = [[]]
        with open(fname, 'r') as file:
            for line in file:
                line = line.split()
                ss.append(line)
        del ss[0:3]
        hierarchy = SspredFile()
        hierarchy.fname = fname
        hierarchy.df = pd.DataFrame(ss, columns=['position', 'aa', 'prediction', '1', '2', '3'])
        return hierarchy

    @staticmethod
    def write(fname):
        raise NotImplementedError('This is not implemented!')


if __name__ == '__main__':

    ss_prediction = PsipredParser().read('/Users/shahrammesdaghi/cmap_plus/w9dy28.ss2')
    for pred in ss_prediction["prediction"]:
        print(pred)
    for aa in ss_prediction["aa"]:
        print(aa)



