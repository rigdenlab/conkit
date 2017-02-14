"""
A module to produce a contact map chord diagram
"""

from __future__ import division

__author__ = "Felix Simkovic"
__date__ = "13 Feb 2017"
__version__ = 0.1

import matplotlib.patches
import matplotlib.pyplot
import numpy

from conkit.plot._Figure import Figure


class ContactMapChordFigure(Figure):
    """A Figure object specifically for a Contact Map chord diagram

    Description
    -----------
    This figure will illustrate the contacts linking the residues
    in the target sequence. This plot is a very common representation
    of contacts. With this figure, you can illustrate intra-molecular.

    Attributes
    ----------
    hierarchy : :obj:`ContactMap <conkit.core.ContactMap>`
       The default contact map hierarchy
    use_conf : bool
       The marker size will correspond to the raw score [default: False]

    Examples
    --------
    >>> import conkit
    >>> cmap = conkit.io.read('toxd/toxd.mat', 'ccmpred').top_map
    >>> conkit.plot.ContactMapChordFigure(cmap)

    """
    AA_ENCODING = {
        'A': '#882D17', 'C': '#F3C300', 'D': '#875692', 'E': '#F38400',
        'F': '#A1CAF1', 'G': '#BE0032', 'H': '#C2B280', 'I': '#848482',
        'K': '#008856', 'L': '#E68FAC', 'M': '#0067A5', 'N': '#F99379',
        'P': '#604E97', 'Q': '#F6A600', 'R': '#B3446C', 'S': '#DCD300',
        'T': '#8DB600', 'V': '#654522', 'W': '#E25822', 'Y': '#2B3D26',
        'X': '#000000'
    }

    def __init__(self, hierarchy, use_conf=False, **kwargs):
        """A new contact map plot

        Parameters
        ----------
        hierarchy : :obj:`ContactMap <conkit.core.ContactMap>`
           The default contact map hierarchy
        use_conf : bool, optional
           The marker size will correspond to the raw score [default: False]
        **kwargs
           General :obj:`Figure <conkit.plot._Figure.Figure>` keyword arguments

        """
        super(ContactMapChordFigure, self).__init__(**kwargs)

        self._hierarchy = None
        self.hierarchy = hierarchy

        self.use_conf = use_conf

        self._draw()

    def __repr__(self):
        return "ContactMapChordFigure(file_name=\"{0}\")".format(self.file_name)

    @property
    def hierarchy(self):
        """The default contact map hierarchy"""
        return self._hierarchy

    @hierarchy.setter
    def hierarchy(self, hierarchy):
        """Define the default contact map hierarchy"""
        if hierarchy:
            Figure._check_hierarchy(hierarchy, "ContactMap")
        if hierarchy.sequence:
            Figure._check_hierarchy(hierarchy.sequence, "Sequence")
        self._hierarchy = hierarchy

    def redraw(self):
        """Re-draw the plot with updated parameters"""
        self._draw()

    def _draw(self):
        """Draw the actual plot"""

        # Re-normalize the data for the lines
        hierarchy = self.hierarchy.rescale()

        # Obtain the data from the hierarchy
        self_data = numpy.asarray([(c.res1, c.res1_seq, c.res2, c.res2_seq, c.raw_score)
                                   for c in hierarchy])
        min_max_data = numpy.append(self_data.T[1], self_data.T[3],).astype(numpy.int64)
        self_data_range = numpy.arange(min_max_data.min(), min_max_data.max() + 1)

        # The number of points in the outer circle
        npoints = self_data_range.shape[0]

        # Calculate the xy coords for each point on the circle
        space = 2 * numpy.pi / npoints
        verts = numpy.zeros((npoints, 2))
        for i in numpy.arange(npoints):
            verts[i] = [npoints * numpy.cos(space * i), npoints * numpy.sin(space * i)]

        # Instantiate the figure
        fig, ax = matplotlib.pyplot.subplots()

        # Calculate and plot the Bezier curves
        bezier_path = numpy.arange(0, 1.01, 0.01)
        for c in self_data:
            x1, y1 = verts[int(c[1]) - self_data_range.min()]
            x2, y2 = verts[int(c[3]) - self_data_range.min()]
            xb, yb = [0, 0]                  # Midpoint the curve is supposed to approach
            x = (1 - bezier_path) ** 2 * x1 + 2 * (1 - bezier_path) * bezier_path * xb + bezier_path ** 2 * x2
            y = (1 - bezier_path) ** 2 * y1 + 2 * (1 - bezier_path) * bezier_path * yb + bezier_path ** 2 * y2
            ax.plot(x, y, color="#000000", alpha=float(c[4]), linestyle="-", zorder=0)

        # Get the amino acids if available
        #     - get the residue data from the original data array
        residue_data = numpy.append(self_data[:, [1, 0]], self_data[:, [3, 2]]).reshape(self_data.T[0].shape[0] * 2, 2)

        #     - compute a default color list
        color_codes = {
            k: ContactMapChordFigure.AA_ENCODING['X']
            for k in self_data_range
        }

        #     - fill default list with data we have
        for k, v in numpy.vstack({tuple(row) for row in residue_data}):
            color_codes[int(k)] = ContactMapChordFigure.AA_ENCODING[v]

        #     - create a color list
        colors = [color_codes[k] for k in sorted(color_codes.keys())]

        # Plot the residue points
        ax.scatter(verts.T[0], verts.T[1], marker='o', color=colors, edgecolors="#000000", zorder=1)

        # # Annotate every residue
        # for c in numpy.vstack({tuple(row) for row in residue_data}):
        #     i = int(c[0]) - self_data_range.min()
        #     angle = numpy.rad2deg(space) * i
        #     xy = verts[i]
        #     label = "{0} {1}".format(c[0], c[1])
        #     ax.annotate(label, xy=xy, xytext=xy, rotation=angle)

        # Prettify the plot
        ax.set_xlim(-npoints - 10, npoints + 11)
        ax.set_ylim(-npoints - 10, npoints + 11)
        ax.axis("off")

        # Make both axes identical in length and remove whitespace around the plot
        aspectratio = Figure._correct_aspect(ax, 1.0)
        ax.set(aspect=aspectratio)
        fig.tight_layout()

        fig.savefig(self.file_name, bbox_inches='tight', dpi=self.dpi)
