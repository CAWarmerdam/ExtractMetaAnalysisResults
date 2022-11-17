#!/usr/bin/env python3

"""
Created:      13/07/2022
Author:       C.A. (Robert) Warmerdam

Copyright (C) 2022 C.A. Warmerdam

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
import sys
import os
import glob

# Metadata
__program__ = "Generate SQLite commands for generating vtables"
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


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    # List all parquet files in the given directory
    virtual_table_creation_format = "CREATE VIRTUAL TABLE eqtls_{} USING parquet('{}')"
    select_virtual_table_format = "SELECT * FROM eqtls_{}"

    virtual_table_creation_commands = list()
    virtual_table_select_commands = list()

    for file_number, parquet_file in enumerate(glob.glob("{}*.parquet".format(argv[0]))):
        virtual_table_creation_commands.append(
            virtual_table_creation_format.format(file_number, parquet_file))

        virtual_table_select_commands.append(
            select_virtual_table_format.format(file_number))

    virtual_table_commands = os.linesep.join(
        virtual_table_creation_commands)
    view_generation_command = os.linesep.join([
        "CREATE VIEW eqtls AS ",
        "{}UNION ".format(os.linesep).join(virtual_table_select_commands)
        ])

    print(virtual_table_commands)
    print(view_generation_command)

    return 0


if __name__ == "__main__":
    sys.exit(main())
