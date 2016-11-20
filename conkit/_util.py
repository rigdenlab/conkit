"""Internal utility functions"""

__author__ = "Felix Simkovic"
__date__ = "20 Nov 2016"
__version__ = "0.1"

import tempfile

def create_tmp_f(mode='w', content=None):
    """Create a temporary file

    Parameters
    ----------
    mode : str, optional
       Write-mode [default: 'w']
    content : str, optional
       The content to write to the file

       If `None`, leave the file empty

    Returns
    -------
    str
       The path to the filename

    """
    f_in = tempfile.NamedTemporaryFile(mode=mode, delete=False)
    if content:
        f_in.write(content)
    f_in.close()
    return f_in.name

