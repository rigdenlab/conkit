from conkit.io._parser import PredictionFileParser
from conkit.core.predictionfile import PredictionFile
from conkit.core.predres import ResiduePrediction
import pandas as pd

class PsipredParser(PredictionFileParser):

    """Class for parsing psipred files."""

    def read(self, f_handle, f_id="psipred"):
        """Return psipred prediction instance."""

        hierarchy = PredictionFile(f_id)

        ss = [[]]

        for line in f_handle:
            line = line.split()
            ss.append(line)
        del ss[0:3]

        df = pd.DataFrame(ss, columns=['position', 'aa', 'prediction', '1', '2', '3'])

        ss_pred = df['prediction'].tolist()


        for idx, x in enumerate(ss_pred):
            residue = ResiduePrediction(str(idx + 1))
            residue.ss2_prediction = x

            hierarchy.add(residue)

        return hierarchy


    @staticmethod
    def write(fname):
        raise NotImplementedError('This is not implemented!')


if __name__ == '__main__':

    import conkit.io

    ss_prediction = conkit.io.read('/Users/shahrammesdaghi/Downloads/w9dy28.ss2', "psipred")
    for residue in ss_prediction:
        print(residue.res_seq)
        print(residue.ss2_prediction)
        print(residue.id)


    # res_name = conkit.io.read('/Users/shahrammesdaghi/Downloads/w9dy28.ss2', "psipred")
    # for residue in res_name:
    #     print(residue.residue_name)




