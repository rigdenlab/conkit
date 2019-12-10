from conkit.io._parser import PredictionFileParser
from conkit.core.predictionfile import PredictionFile
from conkit.core.residueprediction import ResiduePrediction

class iupred2aParser(PredictionFileParser):
    """Class for parsing iupred2a files."""

    def read(self, f_handle, f_id="topcons"):
        """Return iupred2a prediction instance."""

        hierarchy = PredictionFile(f_id)

        for line in f_handle:
            line = line.split()
            if len(line) >= 1:
                residue_position = line[0]
                if residue_position.isnumeric():
                    disorder_score = line[2]
                    residue = ResiduePrediction(residue_position)
                    residue.disorder_prediction = float(disorder_score)
                    hierarchy.add(residue)

        return hierarchy

    @staticmethod
    def write(self):
        raise NotImplementedError('Conkit does not support iupred2a writing!')


