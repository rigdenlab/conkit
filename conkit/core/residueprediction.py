from conkit.core.entity import Entity


class ResiduePrediction(Entity):
    __slots__ = [
        "prediction",
        "res_seq"
    ]

    def __init__(self, res_seq):
        self.res_seq = res_seq
        self.prediction = None

        super(ResiduePrediction, self).__init__(res_seq)
