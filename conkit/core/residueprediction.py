from conkit.core.entity import Entity


class ResiduePrediction(Entity):

    __slots__ = [
        "_prediction",
        "res_seq"
    ]

    def __init__(self, res_seq):

        self.res_seq = res_seq
        self._prediction=None

        super(ResiduePrediction, self).__init__(res_seq)

    @property
    def prediction(self):
        return self._prediction

    @prediction.setter
    def prediction(self, value):
        self._prediction = value

