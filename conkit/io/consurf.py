from conkit.io._parser import PredictionFileParser
from conkit.core.predictionfile import PredictionFile
from conkit.core.residueprediction import ResiduePrediction


class ConsurfParser(PredictionFileParser):
    """Class for parsing consurf files."""

    def read(self, f_handle, f_id="consurf"):
        """Return consurf conservation scores instance."""

        hierarchy = PredictionFile(f_id, predtype = 'conservationpred')

        for line in f_handle:
            line = line.split()
            if len(line) >= 1 and line[0].isnumeric():
                residue = ResiduePrediction(line[0])
                try:
                    residue.prediction = int(line[3].replace('*', ''))
                    hierarchy.add(residue)
                except:
                    print('Check your consurf file')
                    exit()

        return hierarchy

    def write(self):

        raise NotImplementedError('Conkit does not support consurf writing!')


