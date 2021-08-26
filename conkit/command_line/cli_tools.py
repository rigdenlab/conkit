import os
import shutil


def check_file_exists(input_path):
    """Check if a given path exists

    Parameters
    ----------
    input_path : str, None
       Location of the file to be tested

    Returns
    -------
    abspath : str, None
       The absolute path of the file if it exists, None if the input is None

    Raises
    ------
    :exc:`FileNotFoundError`
        The file doesn't exist
    """

    if input_path is None:
        return None
    if os.path.isfile(os.path.abspath(input_path)):
        return os.path.abspath(input_path)
    else:
        raise FileNotFoundError("{} cannot be found".format(input_path))


def prepare_output_directory(outdir, overwrite=False):
    """Prepare the output directory for conkit-validate.

    Parameters
    ----------
    outdir : str
       Path to the output directory
    overwrite : bool
       Whether the output directory should be overwritten or not [default: False]

    Raises
    ------
    :exc:`ValueError`
        The output directory already exists and overwrite is False
    """

    if os.path.isdir(outdir):
        if not overwrite:
            raise ValueError('Output directory already exists')
        else:
            shutil.rmtree(outdir)
    os.mkdir(outdir)


def parse_map_align_stdout(stdout):
    """Parse the stdout of map_align and extract the alignment of residues.

    Parameters
    ----------
    stdout : str
       Starndard output created with map_align

    Returns
    ------
    alignment: dict
        A dictionary where the aligned residue numbers of map_b are the keys and the residue numbers of map_a the values
    """

    alignment = {}
    for line in stdout.split('\n'):
        if line and line.split()[0] == "MAX":
            line = line.rstrip().lstrip().split()
            for residue in line[8:]:
                if residue.count(':') == 1:
                    resnum = residue.split(":")
                    alignment[int(resnum[1])] = int(resnum[0])
                else:
                    print(residue)


    return alignment
