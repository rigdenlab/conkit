# coding=utf-8
#
# BSD 3-Clause License
#
# Copyright (c) 2016-21, University of Liverpool
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
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
"""
Command line object for map_align contact map alignment application
"""

import subprocess


class MapAlignCommandline:
    """
    Command line object for map_align.

    https://github.com/sokrypton/map_align

    Examples
    --------
    >>> cline = MapAlignCommandline(contact_map_a='a.mapalign', contact_map_b='b.mapalign')
    >>> stdout, stderr = cline()

    """

    def __init__(self, cmd="map_align", contact_map_a=None, contact_map_b=None,
                 gap_opening_penalty=None, gap_extension_penalty=None,
                 seq_separation_cutoff=None, n_iterations=None):
        self.cmd = cmd
        self.contact_map_a = contact_map_a
        self.contact_map_b = contact_map_b
        self.gap_opening_penalty = gap_opening_penalty
        self.gap_extension_penalty = gap_extension_penalty
        self.seq_separation_cutoff = seq_separation_cutoff
        self.n_iterations = n_iterations

    def _build_command(self):
        args = [self.cmd, "-a", self.contact_map_a, "-b", self.contact_map_b]
        if self.gap_opening_penalty is not None:
            args += ["-gap_o", str(self.gap_opening_penalty)]
        if self.gap_extension_penalty is not None:
            args += ["-gap_e", str(self.gap_extension_penalty)]
        if self.seq_separation_cutoff is not None:
            args += ["-sep_cut", str(self.seq_separation_cutoff)]
        if self.n_iterations is not None:
            args += ["-iter", str(self.n_iterations)]
        return args

    def __call__(self):
        result = subprocess.run(
            self._build_command(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"map_align exited with code {result.returncode}.\n"
                f"stderr: {result.stderr.strip()}"
            )
        return result.stdout, result.stderr

    def __str__(self):
        return " ".join(self._build_command())
