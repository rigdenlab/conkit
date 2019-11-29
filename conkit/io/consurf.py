from conkit.io._parser import PredictionFileParser
from conkit.core.predictionfile import PredictionFile
from conkit.core.predres import ResiduePrediction


class ConsurfParser(PredictionFileParser):
    """Class for parsing consurf files."""

    def read(self, f_handle, f_id="consurf"):
        """Return consurf conservation scores instance."""

        hierarchy = PredictionFile(f_id)

        tmp_list = []

        # TODO This could be done in one loop

        for line in f_handle:
                line = line.split()
                try:
                    test = line[0]
                    if test.isnumeric():
                        tmp_list.append(line)
                # TODO Catch exeption with if instead
                except:
                    continue

        for idx, x in enumerate(tmp_list):
            # allows both types of consurf output to be parsed
            try:
                x = x[3].replace('*', '')
                x = int(x)
            except:
                x = x[4].replace('*', '')
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