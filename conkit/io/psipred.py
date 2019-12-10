from conkit.io._parser import PredictionFileParser
from conkit.core.predictionfile import PredictionFile
from conkit.core.residueprediction import ResiduePrediction

class PsipredParser(PredictionFileParser):

    """Class for parsing psipred files."""

    def read(self, f_handle, f_id="psipred"):
        """Return psipred prediction instance."""

        hierarchy = PredictionFile(f_id)

        for line in f_handle:
            line = line.split()
            if len(line) >= 1:
                residue_position = line[0]
                if residue_position.isnumeric():
                    ss2_prediction = line[2]
                    residue = ResiduePrediction(residue_position)
                    residue.ss2_prediction = ss2_prediction
                    hierarchy.add(residue)

        return hierarchy


    @staticmethod
    def write(fname):
        raise NotImplementedError('Conkit does not support psipred writing!')





