from conkit.io._parser import PredictionFileParser
from conkit.core.predictionfile import PredictionFile
from conkit.core.residueprediction import ResiduePrediction


class TopconsParser(PredictionFileParser):
    """Class for parsing topcons files."""


    def read(self, f_handle, f_id="topcons"):
        """Return topcons prediction instance."""

        hierarchy = PredictionFile(f_id, predtype='mempred')

        lines = f_handle.readlines()
        membrane_pred = lines[lines.index('TOPCONS predicted topology:\n') + 1]
        membrane_pred_list = list(membrane_pred.rstrip())

        for idx, x in enumerate(membrane_pred_list):
            residue = ResiduePrediction(str(idx + 1))
            residue.prediction = x
            hierarchy.add(residue)

        return hierarchy


    @staticmethod
    def write(self):
        raise NotImplementedError('Conkit does not support topcons writing!')






