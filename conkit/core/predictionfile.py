from conkit.core.entity import Entity


class PredictionFile(Entity):
    """Class to contain residue prediction data"""

    def __init__(self, id, predtype):

        self._predtype=predtype
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


    def add(self, entity):
        """Overwrites conkit.core.entity add method"""

        if not self._is_allowed_prediction_symbol(entity.prediction):
            raise ValueError('Residue prediction symbol {} is not supported for {}'.format(entity.prediction,
                                                                                           self.predtype))
        super(PredictionFile, self).add(entity)



    def _is_allowed_prediction_symbol(self, symbol):

        if self.predtype == 'mempred':
            if symbol not in ['i', 'M', 'o']:
                return False
            else:
                return True
        elif self.predtype == 'sspred':
            if symbol not in ['C', 'E', 'H']:
                return False
            else:
                return True
        elif self.predtype == 'conservationpred':
            if isinstance(symbol, int):
                return True
            else:
                return False
        elif self.predtype == 'dispred':
            if isinstance(symbol, float):
                return True
            else:
                return False


    def search(self, value=None):
        if len(self) == 0:
            return None

        rslt=[]
        if self.predtype in ['mempred', 'sspred','conservationpred']:
            rslt = [residue for residue in self if residue.prediction==value]
        elif self.predtype in ['dispred']:
            rslt = [residue for residue in self if residue.prediction >= value]
        return rslt
