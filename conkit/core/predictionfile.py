from conkit.core.entity import Entity


class PredictionFile(Entity):
    """Class to contain residue prediction data"""

    def __init__(self, id, predtype):

        self._predtype = predtype
        super(PredictionFile, self).__init__(id)

    def __repr__(self):
        return '{}(id="{}" nres={})'.format(self.__class__.__name__, self.id, len(self))

    @property
    def predtype(self):
        return self._predtype

    @predtype.setter
    def predtype(self, value):

        if value not in ['mempred', 'sspred', 'conservationpred', 'dispred']:
            raise ValueError('Supported prediction type are mempred, sspred, conservpred or dispred')

        self._predtype = value

    @property
    def allowed_symbols(self):

        return {'mempred': ['i', 'M', 'o'],
                'sspred': ['C', 'E', 'H'],
                'conservationpred': int,
                'dispred': float
                }

    def add(self, entity):
        """Overwrites conkit.core.entity add method"""

        if not self._is_allowed_prediction_symbol(entity.prediction):
            raise ValueError('Residue prediction symbol {} is not supported for {}'.format(entity.prediction,
                                                                                           self.predtype))
        super(PredictionFile, self).add(entity)

    def _is_allowed_prediction_symbol(self, symbol):
        if isinstance(self.allowed_symbols[self.predtype], list):
            if symbol not in self.allowed_symbols[self.predtype]:
                return False
            else:
                return True
        else:
            if isinstance(symbol, self.allowed_symbols[self.predtype]):
                return True
            else:
                return False

    def search(self, value=None):
        if len(self) == 0:
            return None

        if self.predtype in ['mempred', 'sspred', 'conservationpred']:
            rslt = [residue for residue in self if residue.prediction == value]
            return rslt
        elif self.predtype == 'dispred':
            rslt = [residue for residue in self if residue.prediction >= value]
            return rslt
        else:
            raise ValueError('Unsupported format {}'.format(self.predtype))
