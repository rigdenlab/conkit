"""
A module to produce a domain boundary plot
"""

from __future__ import division

__author__ = "Felix Simkovic"
__date__ = "23 Feb 2017"
__version__ = 0.1

import matplotlib.pyplot
import numpy
import warnings

try:
    import scipy.signal
    SCIPY = True
except ImportError:
    SCIPY = False

try:
    import sklearn.neighbors
    SKLEARN = True
except ImportError:
    SKLEARN = False



from conkit.plot._Figure import Figure
from conkit.plot._plottools import ColorDefinitions


class ContactDensityFigure(Figure):
    """A Figure object specifically for a contact density illustration.

    This figure is an adaptation of the algorithm published by Sadowski
    (2013) [#]_.

    .. [#] Sadowski M. (2013). Prediction of protein domain boundaries
       from inverse covariances. Proteins 81(2), 253-260.

    Attributes
    ----------
    hierarchy : :obj:`ContactMap <conkit.core.ContactMap>`
       The default contact map hierarchy
    bw_method : str
       The method to estimate the bandwidth [default: bowman]

    Examples
    --------
    >>> import conkit
    >>> cmap = conkit.io.read('toxd/toxd.mat', 'ccmpred').top_map
    >>> conkit.plot.ContactDensityFigure(cmap)

    """
    BW_METHODS = ['bowman']

    def __init__(self, hierarchy, bw_method='bowman', **kwargs):
        """A new contact density plot

        Parameters
        ----------
        hierarchy : :obj:`ContactMap <conkit.core.ContactMap>`
           The default contact map hierarchy
        bw_method : str, optional
           The method to estimate the bandwidth [default: bowman]
        **kwargs
           General :obj:`Figure <conkit.plot._Figure.Figure>` keyword arguments

        """
        super(ContactDensityFigure, self).__init__(**kwargs)
        self._bw_method = None
        self._hierarchy = None

        self.bw_method = bw_method
        self.hierarchy = hierarchy

        self._draw()

    def __repr__(self):
        return "{0}(file_name=\"{1}\" bw_method=\"{2}\")".format(
            self.__class__.__name__, self.file_name, self.bw_method
        )

    @property
    def bw_method(self):
        """The method to estimate the bandwidth"""
        return self._bw_method

    @bw_method.setter
    def bw_method(self, bw_method):
        """Define the method to estimate the bandwidth

        Raises
        ------
        ValueError
           Method not yet defined

        """
        if bw_method not in ContactDensityFigure.BW_METHODS:
            msg = "Bandwidth method not yet implemented: {0}".format(bw_method)
            raise ValueError(msg)
        self._bw_method = bw_method

    @property
    def hierarchy(self):
        """A ConKit :obj:`ContactMap <conkit.core.ContactMap>`"""
        return self._hierarchy

    @hierarchy.setter
    def hierarchy(self, hierarchy):
        """Define the ConKit :obj:`ContactMap <conkit.core.ContactMap>`

        Raises
        ------
        RuntimeError
           The hierarchy is not an contact map

        """
        if hierarchy:
            Figure._check_hierarchy(hierarchy, "ContactMap")
        self._hierarchy = hierarchy

    def redraw(self):
        """Re-draw the plot with updated parameters"""
        self._draw()

    def _draw(self):
        """Draw the actual plot
        
        Raises
        ------
        RuntimeError
           Cannot find SciKit package

        """

        if not SKLEARN:
            raise RuntimeError('Cannot find SciKit package')

        # Compute the relevant data we need
        X = numpy.asarray([i for c in self._hierarchy for i in numpy.arange(c.res1_seq, c.res2_seq)])[:, numpy.newaxis]
        X_plot = numpy.linspace(X.min(), X.max(), X.max() - X.min())[:, numpy.newaxis]

        # Obtain the bandwidth as defined by user method
        if self.bw_method == "bowman":
            bandwidth = ContactDensityFigure.bowman_bandwidth(X)
        else:
            bandwidth = ContactDensityFigure.bowman_bandwidth(X)

        # Estimate the Kernel Density using original data and fit random sample
        kde = sklearn.neighbors.KernelDensity(bandwidth=bandwidth).fit(X)
        dens = numpy.exp(kde.score_samples(X_plot))
        
        # Plot the data
        fig, ax = matplotlib.pyplot.subplots()

        ax.plot(X_plot[:, 0], dens, linestyle="solid",
                color=ColorDefinitions.GENERAL, label="Kernel Density Estimate")

        # Find all local minima
        if SCIPY:
            local_minima_idx = scipy.signal.argrelmin(dens)[0]
            ax.scatter(X_plot[local_minima_idx], dens[local_minima_idx], marker="p",
                       color=ColorDefinitions.MISMATCH, label="Local Minimum")
        else:
            warnings.warn("SciPy not installed - cannot determine local minima")

        # Prettify the plot
        ax.set_xlim(X.min(), X.max())
        ax.set_ylim(0., dens.max())

        ax.set_xlabel('Residue number')
        ax.set_ylabel('Kernel Density Estimate')
        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.,
                  scatterpoints=1)

        # Make axes length proportional and remove whitespace around the plot
        aspectratio = Figure._correct_aspect(ax, 0.3)
        ax.set(aspect=aspectratio)
        fig.tight_layout()

        fig.savefig(self.file_name, bbox_inches='tight', dpi=self.dpi)

    @staticmethod
    def bowman_bandwidth(X):
        """This is the optimal bandwidth if the point distribution is Gaussian.

        To calculate the bandwidth for the 1D data array ``X`` with ``n`` data points, the following
        equation is used:
        
        .. math::

           bandwidth=\\sqrt{\\frac{\\sum{X}^2}{n}-(\\frac{\\sum{X}}{n})^2}*(\\frac{3*n}{4})^\\frac{-1}{5}
        
        This equation is a direct implementation taken from Bowman & Azzalini [#]_.
        
        .. [#] Bowman, A.W. & Azzalini, A. (1997). Applied Smoothing Techniques for Data Analysis.

        Parameters
        ----------
        X : list, tuple
           A list of data points


        """
        X = numpy.asarray(X)
        sigma = numpy.sqrt((X ** 2).sum() / X.shape[0] - (X.sum() / X.shape[0]) ** 2)
        return sigma * ((3. * X.shape[0] / 4.) ** (-1. / 5.))
