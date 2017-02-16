"""Internal utility functions"""

__author__ = "Felix Simkovic"
__date__ = "16 Feb 2017"
__version__ = 0.1

import numpy


def points_on_circle(radius, h=0, k=0):
    """Calculate points on a circle with even spacing

    Parameters
    ----------
    radius : int
       The radius of the circle
    h : int, optional
       The x coordinate of the origin
    k : int, optional
       The y coordinate of the origin

    Returns
    -------
    list
       The list of coordinates for each point

    """
    space = 2 * numpy.pi / radius
    coords = numpy.zeros((radius, 2))
    for i in numpy.arange(radius):
        coords[i] = [
            round(h + radius * numpy.cos(space * i), 6),
            round(k + radius * numpy.sin(space * i), 6)
        ]
    return coords.tolist()
