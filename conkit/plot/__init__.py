"""Plot interface for automated figure generation"""

__author__ = "Felix Simkovic"
__date__ = "07 Feb 2017"
__version__ = 0.1

import matplotlib
matplotlib.use('Agg')

from conkit.plot.ContactMapPlot import ContactMapFigure
from conkit.plot.ContactMapChordPlot import ContactMapChordFigure
from conkit.plot.PrecisionEvaluationPlot import PrecisionEvaluationFigure
from conkit.plot.SequenceCoveragePlot import SequenceCoverageFigure
