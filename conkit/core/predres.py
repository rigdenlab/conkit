from conkit.core.entity import Entity
import numpy as np


class ResiduePrediction(Entity):

    __slots__ = [
        "_conservation_score",
        "_ss2_prediction",
        "_membrane_prediction",
        "_disorder_prediction",
        "res_seq"
    ]

    def __init__(self, res_seq):

        self.res_seq = res_seq

        self._conservation_score = np.nan
        self._ss2_prediction = np.nan
        self._membrane_prediction = np.nan
        self._disorder_prediction = np.nan

        super(ResiduePrediction, self).__init__(res_seq)


    @property
    def ss2_prediction(self):
        return self._ss2_prediction

    @ss2_prediction.setter
    def ss2_prediction(self, value):
        self._ss2_prediction = value

    @property
    def conservation_score(self):
        return self._conservation_score

    @conservation_score.setter
    def conservation_score(self, value):
        self._conservation_score = value

    @property
    def membrane_prediction(self):
        return self._membrane_prediction

    @membrane_prediction.setter
    def membrane_prediction(self, value):
        self._membrane_prediction = value

    @property
    def disorder_prediction(self):
        return self._disorder_prediction

    @disorder_prediction.setter
    def disorder_prediction(self, value):
        self._disorder_prediction = value

