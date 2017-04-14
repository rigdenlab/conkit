"""Command line facility for ConKit scripts"""

__author__ = "Felix Simkovic"
__date__ = "14 Apr 2017"
__version__ = "0.1"

import logging


def get_logger(name, level='info'):
    """Return an initialized logging handler

    Parameters
    ----------
    name : str
       The logger name
    level : str, optional
       The logging level to be used [default: info]
       To change, use one of 
           [ notset | info | debug | warning | error | critical ]
    
    Returns
    -------
    logger
       Instance of a :obj:`logger <logging.Logger>`

    Raises
    ------
    ValueError
       Unknown logging level

    """
    # Define the logging level handlers
    logger_levels = {'notset': logging.NOTSET, 'info': logging.INFO, 'debug': logging.DEBUG,
                     'warning': logging.WARNING, 'error': logging.ERROR, 'critical': logging.CRITICAL}        
    if level in logger_levels:
        LEVEL = logger_levels[level]
    else:
        raise ValueError("Unknown logging level: %s" % level)

    FORMAT = "%(message)s" 
    logging.basicConfig(format=FORMAT, level=LEVEL)
    return logging.getLogger(name)

