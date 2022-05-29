# BSD 3-Clause License
#
# Copyright (c) 2016-21, University of Liverpool
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""Plot interface for automated figure generation"""

__author__ = "Felix Simkovic"
__date__ = "07 Feb 2017"
__version__ = "0.13.3"

import matplotlib

matplotlib.use("Agg")


def ContactMapFigure(*args, **kwargs):
    """:obj:`~conkit.plot.contactmap.ContactMapFigure` instance"""
    from conkit.plot.contactmap import ContactMapFigure

    return ContactMapFigure(*args, **kwargs)


def ContactMapChordFigure(*args, **kwargs):
    """:obj:`~conkit.plot.contactmapchord.ContactMapChordFigure` instance"""
    from conkit.plot.contactmapchord import ContactMapChordFigure

    return ContactMapChordFigure(*args, **kwargs)


def ContactMapMatrixFigure(*args, **kwargs):
    """:obj:`~conkit.plot.contactmatrix.ContactMapMatrixFigure` instance"""
    from conkit.plot.contactmapmatrix import ContactMapMatrixFigure

    return ContactMapMatrixFigure(*args, **kwargs)


def ContactDensityFigure(*args, **kwargs):
    """:obj:`~conkit.plot.contactdensity.ContactDensityFigure` instance"""
    from conkit.plot.contactdensity import ContactDensityFigure

    return ContactDensityFigure(*args, **kwargs)


def PrecisionEvaluationFigure(*args, **kwargs):
    """:obj:`~conkit.plot.precisionevaluation.PrecisionEvaluationFigure` instance"""
    from conkit.plot.precisionevaluation import PrecisionEvaluationFigure

    return PrecisionEvaluationFigure(*args, **kwargs)


def SequenceCoverageFigure(*args, **kwargs):
    """:obj:`~conkit.plot.sequencecoverage.SequenceCoverageFigure` instance"""
    from conkit.plot.sequencecoverage import SequenceCoverageFigure

    return SequenceCoverageFigure(*args, **kwargs)


def DistogramHeatmapFigure(*args, **kwargs):
    """:obj:`~conkit.plot.distogramheatmap.DistogramHeatmapFigure` instance"""
    from conkit.plot.distogramheatmap import DistogramHeatmapFigure

    return DistogramHeatmapFigure(*args, **kwargs)


def ModelValidationFigure(*args, **kwargs):
    """:obj:`~conkit.plot.modelvalidation.ModelValidationFigure` instance"""
    from conkit.plot.modelvalidation import ModelValidationFigure

    return ModelValidationFigure(*args, **kwargs)
