# BSD 3-Clause License
#
# Copyright (c) 2016-18, University of Liverpool
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice, #   this list of conditions and the following disclaimer in the documentation
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
"""Various miscellaneous code required by ConKit"""

__author__ = "Felix Simkovic"
__date__ = "18 May 2018"
__version__ = "2.0"

import numpy as np
import warnings


def fAND(a, b):
    """Function AND operator for two arguments"""
    return a and b


def fOR(a, b):
    """Function OR operator for two arguments"""
    return a or b


def deprecate(version, msg=None):
    """Decorator to deprecate Python classes and functions
    
    Parameters
    ----------
    version : str
       A string containing the version with which the callable is removed
    msg : str, optional
       An additional message that will be displayed alongside the default message

    
    Examples
    --------
    Enable :obj:`~DeprecationWarning` messages to be displayed.

    >>> import warnings
    >>> warnings.simplefilter('default')
    
    Decorate a simple Python function without additional message
    
    >>> @deprecate('0.0.0')
    ... def sum(a, b):
    ...     return a + b
    >>> sum(1, 2)
    deprecated.py:34: DeprecationWarning: sum has been deprecated and will be removed in version 0.0.0!
      warnings.warn(message, DeprecationWarning)
    3

    Decorate a simple Python function with additional message

    >>> @deprecate('0.0.1', msg='Use XXX instead!')
    ... def sum(a, b):
    ...     return a + b
    >>> sum(2, 2)
    deprecated.py:34: DeprecationWarning: sum has been deprecated and will be removed in version 0.0.0! - Use XXX instead!
      warnings.warn(message, DeprecationWarning)
    4

    Decorate an entire Python class

    >>> @deprecate('0.0.2')
    ... class Obj(object):
    ...     pass
    >>> Obj()
    deprecated.py:34: DeprecationWarning: Obj has been deprecated and will be removed in version 0.0.2!
      warnings.warn(message, DeprecationWarning)
    <__main__.Obj object at 0x7f8ee0f1ead0>

    Decorate a Python class method

    >>> class Obj(object):
    ...     def __init__(self, v):
    ...         self.v = v
    ...     @deprecate('0.0.3')
    ...     def mul(self, other):
    ...         return self.v * other.v
    >>> Obj(2).mul(Obj(3))
    deprecated.py:34: DeprecationWarning: mul has been deprecated and will be removed in version 0.0.3!
      warnings.warn(message, DeprecationWarning)
    6

    Decorate a Python class staticmethod
    
    >>> class Obj(object):
    ...     @staticmethod
    ...     @deprecate('0.0.4')
    ...     def sub(a, b):
    ...         return a - b
    ... 
    >>> Obj.sub(2, 1)
    deprecated.py:34: DeprecationWarning: sub has been deprecated and will be removed in version 0.0.4!
      warnings.warn(message, DeprecationWarning)
    1

    Decorate a Python class classmethod

    >>> class Obj(object):
    ...     CONST = 5
    ...     @classmethod
    ...     @deprecate('0.0.5')
    ...     def sub(cls, a):
    ...         return a - cls.CONST
    ... 
    >>> Obj().sub(5)
    deprecated.py:34: DeprecationWarning: sub has been deprecated and will be removed in version 0.0.5!
      warnings.warn(message, DeprecationWarning)
    0

    """

    def deprecate_decorator(callable_):
        def warn(*args, **kwargs):
            message = "%s has been deprecated and will be removed in version %s!" % (callable_.__name__, version)
            if msg:
                message += " - %s" % msg
            warnings.warn(message, DeprecationWarning)
            return callable_(*args, **kwargs)

        return warn

    return deprecate_decorator


def normalize(data, vmin=0, vmax=1):
    """Apply a Feature scaling algorithm to normalize the data

    This normalization will bring all values into the range [0, 1]. This function
    allows range restrictions by values ``vmin`` and ``vmax``.

    .. math::

       {X}'=\\frac{(X-X_{min})(vmax-vmin)}{X_{max}-X_{min}}

    Parameters
    ----------
    data : list, tuple
       The data to normalize
    vmin : int, float, optional
       The minimum value
    vmax : int, float, optional
       The maximum value

    Returns
    -------
    list
       The normalized data

    """
    data = np.array(data, dtype=np.float64)
    return (vmin + (data - data.min()) * (vmax - vmin) / (data.max() - data.min())).tolist()
