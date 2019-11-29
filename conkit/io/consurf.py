from conkit.io._parser import PredictionFileParser
from conkit.core.predictionfile import PredictionFile
from conkit.core.predres import ResiduePrediction


class ConsurfParser(PredictionFileParser):
    """Class for parsing consurf files."""

    def read(self, f_handle, f_id="consurf"):
        """Return consurf conservation scores instance."""

        hierarchy = PredictionFile(f_id)

        tmp_list = []

        for line in f_handle:
                line = line.split()
                try:
                    test = line[0]
                    if test.isnumeric():
                        tmp_list.append(line)
                except Exception as e:
                    if str(e) == 'list index out of range':
                        continue

        for idx, x in enumerate(tmp_list):
            x = x[3].replace('*', '')
            x = int(x)
            residue = ResiduePrediction(str(idx+1))
            residue.conservation_score = x

            hierarchy.add(residue)

        return hierarchy

    def write(self):

        raise NotImplementedError('This is not implemented!')


if __name__ == "__main__":

    import conkit.io

    conservation_prediction = conkit.io.read('/Users/shahrammesdaghi/Downloads/consurf.grades', "consurf")


    for residue in conservation_prediction:
        #print(residue.res_seq)
        print(residue.conservation_score)
        #print(residue.id)