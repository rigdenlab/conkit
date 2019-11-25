import conkit.io
import pandas as pd

class ContactPlotPlus(object):

    def __init__(self, contactmap):
        self._ss_pred = None
        self._con_pred = None

    @property
    def con_pred(self):
        return self._con_pred

    @con_pred.setter
    def con_pred(self, value):
        self._con_pred = value

    @property
    def ss_pred(self):
        return self._ss_pred

    @ss_pred.setter
    def ss_pred(self, value):
        self._ss_pred = value

    @property
    def layers_to_plot(self):
        rslt = []

        for layer in [self.ss_pred, self.con_pred]:
            if layer is not None:
                rslt.append(layer)
        return rslt

    def plot(self):
        consurf_prediction = conkit.io.ConsurfParser().read('/Users/shahrammesdaghi/cmap_plus/consurf.grades')
        consurf = consurf_prediction.con_score

        ss_prediction = conkit.io.PsipredParser().read('/Users/shahrammesdaghi/cmap_plus/w9dy28.ss2')
        ss = ss_prediction.df


        print(self.layers_to_plot)


#my_plot(conkit.map)
#my_plot.ss_pred = SspredFile()
#my_plot.plot()

