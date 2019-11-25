class MempredFile(object):

    """Class to contain TM prediction data from a file"""

    def __init__(self):
        self._fname = None
        self._membrane_pred_list = None

    @property
    def membrane_pred_list(self):
        return self._membrane_pred_list

    @membrane_pred_list.setter
    def membrane_pred_list(self, value):
        self._membrane_pred_list = value

    @property
    def fname(self):
        return self._fname

    @fname.setter
    def fname(self, value):
        self._fname = value


class TopconsParser(object):
    """Class for parsing topcons files."""

    def __init__(self):
        """Initialise an instance of the topconsParser."""
        pass

    @staticmethod
    def read(fname):
        """Return topcons prediction instance."""

        with open(fname, 'r') as file:
            topcons = file.readlines()

        topcons_index = topcons.index('TOPCONS predicted topology:\n')
        membrane_pred = topcons[topcons_index + 1]
        membrane_pred_list = []
        for x in membrane_pred:
            membrane_pred_list.append(x)
        membrane_pred_list.remove('\n')
        hierarchy = MempredFile()
        hierarchy.fname = fname
        hierarchy.membrane_pred_list = membrane_pred_list
        return hierarchy

    @staticmethod
    def write(self):
        raise NotImplementedError('This is not implemented!')


if __name__ == '__main__':

    topcons_prediction = TopconsParser().read('/Users/shahrammesdaghi/cmap_plus/query.result.txt')
    print(topcons_prediction.membrane_pred_list)




