"""
Created:      20/02/2024
Author:       C.A. (Robert) Warmerdam

Copyright (C) 2024 C.A. Warmerdam

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""

# Standard imports.
import os
import sys
import argparse

# Metadata
__program__ = "CNV-caller"
__author__ = "C.A. (Robert) Warmerdam"
__email__ = "c.a.warmerdam@umcg.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


# Constants

# Classes

# Functions
import numpy as np


def calculate_ld(dosages):
    alt_allele_frequencies = dosages.sum(axis=1) / 2 * dosages.shape[1]
    ref_allele_frequencies = 1 - alt_allele_frequencies

def calculate_ld(dosages_a, dosages_b):
    return np.corrcoef(dosages_a, dosages_b)


# Main
def main(argv=None):
    if argv is None:
        argv = sys.argv
    # Process input

    with open("correlations.csv", 'w') as opened:
        print("")
        for i in range(dosages.shape[0]):
            snp_a = dosages.index[i]
            row_a = dosages.iloc[i]
            print("\rTest {} / {}".format(i+1, dosages.shape[0]), end="", flush=True, file=sys.stdout)
            start = i+1
            if i+1 < dosages.shape[0]:
                for j in range(start,dosages.shape[0]):
                    snp_b = dosages.index[j]
                    row_b = dosages.iloc[j]
                    corr = np.corrcoef(row_a, row_b)[1, 0]
                    chars = opened.write("{}\t{}\t{}\n".format(snp_a, snp_b, corr))

    # Perform method
    # Output
    return 0


if __name__ == "__main__":
    sys.exit(main())
