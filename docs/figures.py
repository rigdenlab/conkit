__author__ = "Felix Simkovic"

import logging
import os
import subprocess
import shutil

from conkit.command_line.conkit_plot import main as plot_cmd

TMP_DIR = "_tmp"
DATA_DIR = os.path.join(TMP_DIR, "conkit-examples")
OUT_DIR = os.path.join("examples", "images")

PLOTS = [
    [
        "chord",
        "--confidence",
        "--overwrite",
        "-o",
        os.path.join(OUT_DIR, "toxd_chord_simple.png"),
        os.path.join(DATA_DIR, "toxd", "toxd.fasta"),
        "fasta",
        os.path.join(DATA_DIR, "toxd", "toxd.mat"),
        "ccmpred",
    ],
    [
        "cmap",
        "--overwrite",
        "-o",
        os.path.join(OUT_DIR, "toxd_cmap_simple.png"),
        os.path.join(DATA_DIR, "toxd", "toxd.fasta"),
        "fasta",
        os.path.join(DATA_DIR, "toxd", "toxd.mat"),
        "ccmpred",
    ],
    [
        "cmap",
        "--overwrite",
        "-o",
        os.path.join(OUT_DIR, "toxd_cmap_reference.png"),
        "-p",
        os.path.join(DATA_DIR, "toxd", "toxd.pdb"),
        "-pf",
        "pdb",
        os.path.join(DATA_DIR, "toxd", "toxd.fasta"),
        "fasta",
        os.path.join(DATA_DIR, "toxd", "toxd.mat"),
        "ccmpred",
    ],
    [
        "cmap",
        "--overwrite",
        "-o",
        os.path.join(OUT_DIR, "toxd_cmap_advanced.png"),
        "-e",
        os.path.join(DATA_DIR, "toxd", "toxd.psicov"),
        "-ef",
        "psicov",
        "-p",
        os.path.join(DATA_DIR, "toxd", "toxd.pdb"),
        "-pf",
        "pdb",
        os.path.join(DATA_DIR, "toxd", "toxd.fasta"),
        "fasta",
        os.path.join(DATA_DIR, "toxd", "toxd.mat"),
        "ccmpred",
    ],
    [
        "cmap",
        "--overwrite",
        "--confidence",
        "-o",
        os.path.join(OUT_DIR, "toxd_cmap_confidence.png"),
        "-e",
        os.path.join(DATA_DIR, "toxd", "toxd.psicov"),
        "-ef",
        "psicov",
        "-p",
        os.path.join(DATA_DIR, "toxd", "toxd.pdb"),
        "-pf",
        "pdb",
        os.path.join(DATA_DIR, "toxd", "toxd.fasta"),
        "fasta",
        os.path.join(DATA_DIR, "toxd", "toxd.mat"),
        "ccmpred",
    ],
    [
        "cmat",
        "--overwrite",
        "-o",
        os.path.join(OUT_DIR, "toxd_cmat_simple.png"),
        os.path.join(DATA_DIR, "toxd", "toxd.fasta"),
        "fasta",
        os.path.join(DATA_DIR, "toxd", "toxd.mat"),
        "ccmpred",
    ],
    [
        "cmat",
        "--overwrite",
        "-o",
        os.path.join(OUT_DIR, "toxd_cmat_advanced.png"),
        "-e",
        os.path.join(DATA_DIR, "toxd", "toxd.psicov"),
        "-ef",
        "psicov",
        os.path.join(DATA_DIR, "toxd", "toxd.fasta"),
        "fasta",
        os.path.join(DATA_DIR, "toxd", "toxd.mat"),
        "ccmpred",
    ],
    [
        "scov",
        "--overwrite",
        "-o",
        os.path.join(OUT_DIR, "toxd_scov_plot.png"),
        os.path.join(DATA_DIR, "toxd", "toxd.a3m"),
        "a3m",
    ],
    [
        "peval",
        "--overwrite",
        "-o",
        os.path.join(OUT_DIR, "toxd_peval_plot.png"),
        "-j",
        "0.1",
        "-min",
        "0",
        "-max",
        "2",
        os.path.join(DATA_DIR, "toxd", "toxd.pdb"),
        "pdb",
        os.path.join(DATA_DIR, "toxd", "toxd.fasta"),
        "fasta",
        os.path.join(DATA_DIR, "toxd", "toxd.mat"),
        "ccmpred",
    ],
    [
        "cdens",
        "--overwrite",
        "-o",
        os.path.join(OUT_DIR, "4p9g_cdens_plot.png"),
        os.path.join(DATA_DIR, "4p9g", "4p9g.fasta"),
        "fasta",
        os.path.join(DATA_DIR, "4p9g", "4p9g.mat"),
        "ccmpred",
    ],
]


def prepare():
    for directory in {TMP_DIR, "examples/images"}:
        if os.path.exists(directory):
            shutil.rmtree(directory)
        os.makedirs(directory)

    subprocess.call(["git", "clone", "https://github.com/rigdenlab/conkit-examples.git", DATA_DIR])

def generate_all():
    prepare()
    logging.disable(level=logging.CRITICAL)
    for plot_argv in PLOTS:
        generate(plot_argv)
    logging.disable(level=logging.NOTSET)
    cleanup()


def generate(argv):
    plot_cmd(argv)


def cleanup():
    if os.path.exists(TMP_DIR):
        shutil.rmtree(TMP_DIR)

