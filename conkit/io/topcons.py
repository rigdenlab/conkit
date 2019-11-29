from conkit.io._parser import PredictionFileParser
from conkit.core.predictionfile import PredictionFile
from conkit.core.predres import ResiduePrediction


class TopconsParser(PredictionFileParser):
    """Class for parsing topcons files."""


    def read(self, f_handle, f_id="topcons"):
        """Return topcons prediction instance."""

        hierarchy = PredictionFile(f_id)
        topcons_prediction_file = list(f_handle)
        topcons_index = topcons_prediction_file.index('TOPCONS predicted topology:\n')
        membrane_pred = topcons_prediction_file[topcons_index + 1]
        membrane_pred_list = []

        for x in membrane_pred:
            membrane_pred_list.append(x)
        membrane_pred_list.remove('\n')

        for idx, x in enumerate(membrane_pred_list):
            residue = ResiduePrediction(str(idx + 1))
            residue.membrane_prediction = x

            hierarchy.add(residue)

        return hierarchy


    @staticmethod
    def write(self):
        raise NotImplementedError('This is not implemented!')


if __name__ == '__main__':

    import conkit.io
    mem_prediction = conkit.io.read('/Users/shahrammesdaghi/Downloads/query.result.txt', "topcons")
    for residue in mem_prediction:
        #print(residue.res_seq)
        print(residue.membrane_prediction)
        #print(residue.id)




