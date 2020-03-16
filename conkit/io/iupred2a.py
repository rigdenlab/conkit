from conkit.io._parser import PredictionFileParser
from conkit.core.predictionfile import PredictionFile
from conkit.core.residueprediction import ResiduePrediction

class iupred2aParser(PredictionFileParser):
    """Class for parsing iupred2a files."""

    def read(self, f_handle, f_id="topcons"):
        """Return iupred2a prediction instance."""

        hierarchy = PredictionFile(f_id, predtype='dispred')

        for line in f_handle:
            line = line.split()
            if len(line) >= 1 and line[0].isnumeric():
                residue = ResiduePrediction(line[0])
                residue.prediction = float(line[2])
                hierarchy.add(residue)

        return hierarchy

    @staticmethod
    def write(self):
        raise NotImplementedError('Conkit does not support iupred2a writing!')


