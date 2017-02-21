"""Internal utility functions"""

__author__ = "Felix Simkovic"
__date__ = "16 Feb 2017"
__version__ = 0.1

import numpy


class ColorDefinitions(object):
    """A class storing all color definitions for the various plots
    for fast and easy handling
    """
    
    # Contact map colors
    MATCH = '#0F0B2C'
    MISMATCH = '#DC4869'
    STRUCTURAL = '#D8D6D6'

    # Sequence coverage colors
    L5CUTOFF = '#3F4587'
    L20CUTOFF = '#B5DD2B'

    # Precision evaluation colors
    PRECISION50 = L5CUTOFF
    FACTOR1 = L20CUTOFF

    # Chord plot encoding
    AA_ENCODING = { 
        'A': '#882D17', 'C': '#F3C300', 'D': '#875692', 'E': '#F38400',
        'F': '#A1CAF1', 'G': '#BE0032', 'H': '#C2B280', 'I': '#848482',
        'K': '#008856', 'L': '#E68FAC', 'M': '#0067A5', 'N': '#F99379',
        'P': '#604E97', 'Q': '#F6A600', 'R': '#B3446C', 'S': '#DCD300',
        'T': '#8DB600', 'V': '#654522', 'W': '#E25822', 'Y': '#2B3D26',
        'X': '#000000'
    } 
    
    # General
    GENERAL = '#000000'




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

