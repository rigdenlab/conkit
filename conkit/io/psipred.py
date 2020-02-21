from conkit.io._parser import PredictionFileParser
from conkit.core.predictionfile import PredictionFile
from conkit.core.residueprediction import ResiduePrediction

class PsipredParser(PredictionFileParser):

    """Class for parsing psipred files."""

    def read(self, f_handle, f_id="psipred"):
        """Return psipred prediction instance."""

        hierarchy = PredictionFile(f_id, predtype='sspred')

        for line in f_handle:
            line = line.split()
            if len(line) >= 1 and line[0].isnumeric():
                residue = ResiduePrediction(line[0])
                residue.prediction = line[2]
                hierarchy.add(residue)
        return hierarchy

    @staticmethod
    def write(fname):
        raise NotImplementedError('Conkit does not support psipred writing!')




