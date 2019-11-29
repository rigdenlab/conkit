from conkit.io._parser import PredictionFileParser
from conkit.core.predictionfile import PredictionFile
from conkit.core.predres import ResiduePrediction
import pandas as pd


class iupred2aParser(PredictionFileParser):
    """Class for parsing iupred2a files."""

    def read(self, f_handle, f_id="topcons"):
        """Return iupred2a prediction instance."""

        hierarchy = PredictionFile(f_id)

        iupred_data = pd.read_csv(f_handle, sep='\t', comment='#', skip_blank_lines=True, header=None)
        iupred2a = iupred_data[2].tolist()

        for idx, x in enumerate(iupred2a):
            residue = ResiduePrediction(str(idx + 1))
            residue.disorder_prediction = x

            hierarchy.add(residue)

        return hierarchy

    @staticmethod
    def write(self):
        raise NotImplementedError('This is not implemented!')


if __name__ == '__main__':

    import conkit.io

    disorder_prediction = conkit.io.read('/Users/shahrammesdaghi/Downloads/iupred2a.txt', "iupred2a")
    for residue in disorder_prediction:
        # print(residue.res_seq)
        print(residue.disorder_prediction)
        # print(residue.id)