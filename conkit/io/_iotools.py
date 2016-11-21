"""Internal utility functions"""

__author__ = "Felix Simkovic"
__date__ = "20 Nov 2016"
__version__ = 0.1

import sys
import tempfile

MODE_APPENDIX = 'b' if sys.version_info.major == 3 else ''


def create_tmp_f(content=None, mode='w'):
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


def is_str_like(content):
    """Check if an instance is string-like

    Parameters
    ----------
    content : str
       The content to check

    Returns
    -------
    bool
       True if is str like else False

    Notes
    -----
    Function taken from Numpy, credits to the author

    """
    try:
        content + ''
    except (TypeError, ValueError):
        return False
    return True


def open_f_handle(f_handle, mode):
    """Open a filehandle

    Parameters
    ----------
    f_handle : file_handle, file_name
       A file handle or a file name
    mode : str
       read, write or append

    Returns
    -------
    f_handle
       The opened file handle

    Raises
    ------
    TypeError
       f_handle must be str of filehandle
    ValueError
       Mode needs to be one of: append, read, write

    """
    # Check the mode of opening the file
    if mode not in ['append', 'read', 'write']:
        raise ValueError('Mode needs to be one of: append, read, write')
    mode = mode[0]
    
    try:
        if is_str_like(f_handle):
            # Deal with some Python2/3 compatibility issues
            if sys.version_info.major >= 3:
                f_handle = open(f_handle, mode + 'b')
            else:
                f_handle = open(f_handle, mode)
        elif f_handle.mode == mode + MODE_APPENDIX:
            f_handle = f_handle
        else:
            raise TypeError("f_handle must be str of filehandle")
    except AttributeError:
        raise TypeError("f_handle must be str of filehandle")

    return f_handle

