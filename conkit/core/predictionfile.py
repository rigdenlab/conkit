from conkit.core.entity import Entity



class PredictionFile(Entity):

    """Class to contain residue prediction data"""

    def __init__(self, id):
        super(PredictionFile, self).__init__(id)

    def __repr__(self):
        return '{}(id="{}" nres={})'.format(self.__class__.__name__, self.id, len(self))
