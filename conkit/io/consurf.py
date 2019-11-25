class ConsurfFile(object):

    """Class to contain consurf data from a file"""

    def __init__(self):
        self._fname = None
        self._con_score = None

    @property
    def con_score(self):
        return self._con_score

    @con_score.setter
    def con_score(self, value):
        self._con_score = value

    @property
    def fname(self):
        return self._fname

    @fname.setter
    def fname(self, value):
        self._fname = value


class ConsurfParser(object):
    """Class for parsing consurf files."""

    def __init__(self):

        """Initialise an instance of the psipredParser."""

        pass

    @staticmethod
    def read(fname):
        """Return consurf conservation scores instance."""

        consurf = []
        with open(fname, 'r') as file:
            for line in file:
                line = line.split()
                try:
                    test = line[0]
                    if test.isnumeric():
                        consurf.append(line)
                except:
                    continue

        con_score = []
        for x in consurf:
            # allows both types of consurf output to be parsed
            try:
                x = x[3].replace('*','')
                x = int(x)
            except:
                x = x[4].replace('*', '')
                x = int(x)
            con_score.append(x)


        hierarchy = ConsurfFile()
        hierarchy.fname = fname
        hierarchy.con_score = con_score
        return hierarchy

    @staticmethod
    def write(self):
        raise NotImplementedError('This is not implemented!')


if __name__ == '__main__':

    consurf_prediction = ConsurfParser().read('/Users/shahrammesdaghi/cmap_plus/consurf.grades')
    print(consurf_prediction.con_score)
