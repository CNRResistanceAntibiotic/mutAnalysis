#!/usr/bin/env python3
"""
The version is stored here in a separate file so it can exist in only one place.
https://stackoverflow.com/a/7071358/2438989
Copyright 2020 Aurélien BIRER (abirer36@gmail.com)
https://github.com/CNRResistanceAntibiotic/mutAnalysis

This file is part of CGST. CGST is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. CGST is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with CGST. If
not, see <http://www.gnu.org/licenses/>.
"""
import argparse
import os
import re
import shutil

from mutAnalysis import mapping, mut2report
from mutAnalysis import bam2count, utils


def main(args):

    print("Version mutAnalysis: ", version())

    reads_1 = os.path.abspath(args.reads_1)
    reads_2 = os.path.abspath(args.reads_2)
    wk_dir = os.path.abspath(args.workDir)
    initial = args.initial
    force = args.force

    if not reads_1:
        print("Reads R1 file is missing !\n", flush=True)
        print(usage, flush=True)
        exit(1)

    if not reads_2:
        print("Reads R1 file is missing !\n", flush=True)
        print(usage, flush=True)
        exit(1)

    if not wk_dir:
        print("Working directory is missing !\n", flush=True)
        print(usage, flush=True)
        exit(1)

    if not initial:
        print("Initial is missing !\n", flush=True)
        print(usage, flush=True)
        exit(1)

    dir_path = os.path.dirname(os.path.realpath(__file__))

    database = os.path.join(os.path.dirname(dir_path), "database")
    mutation_database = os.path.join(database, "mutations.tsv")
    sequence_file = os.path.join(database, "sequences.fasta")

    if not os.path.exists(wk_dir):
        os.makedirs(wk_dir)

    # print folders/files path
    print("Reads R1 File: {0}".format(reads_1), flush=True)
    print("Reads R2 File: {0}".format(reads_2), flush=True)
    print("Work directory: {0}".format(wk_dir), flush=True)
    print("Initial user: {0}".format(initial), flush=True)
    print("Force: {0}".format(force), flush=True)
    print("Application run at : {0}\n".format(dir_path), flush=True)

    #########################################

    mut_dict = utils.read_mutation_database(mutation_database)

    ##################

    sequence_new_tmp_file = os.path.join(wk_dir, "sequence_tmp.fasta")
    sequence_new_file = os.path.join(wk_dir, "sequence.fasta")
    shutil.copy(sequence_file, sequence_new_tmp_file)
    # rename
    with open(sequence_new_tmp_file, "r") as file:
        with open(sequence_new_file, "w") as out_file:
            for line in file.readlines():
                if ">" in line:
                    line = line.replace(':', '_').replace('(', '').replace(')', '').replace('\'', 'pr').replace('-', '')
                    out_file.write(line)
                else:
                    out_file.write(line)
    os.remove(sequence_new_tmp_file)
    sequence_file = sequence_new_file


    #######################################

    for feature_name, mutation_dict in mut_dict.items():

        print("FEATURE : {0}".format(feature_name), flush=True)

        feature_name = feature_name.replace(':', '_').replace('(', '').replace(')', '').replace('\'', 'pr').replace('-', '')

        print("FEATURE corrected : {0}".format(feature_name), flush=True)

        print("\n-----------------", flush=True)
        print("MAPPING READS ON SEQUENCES", flush=True)
        print("-----------------", flush=True)
        # Mapping
        mapping.main(sequence_file, reads_1, reads_2, wk_dir)

        print("\n-----------------", flush=True)
        print("COUNT MUTATION ON ALIGNEMENT", flush=True)
        print("-----------------", flush=True)

        for mut_prot in mutation_dict["proteic"]:
            pattern = re.compile('([a-zA-Z_-]+)*([0-9]*)([a-zA-Z_-]+)')
            match = pattern.match(mut_prot)
            if match:
                pos_mutation = int(match.groups()[1])
            # Count
            position = "{0}:{1}-{2}".format(feature_name, (pos_mutation*3)-2, (pos_mutation*3))
            bam2count.main(wk_dir, sequence_file, position, feature_name)

            print("\n-----------------", flush=True)
            print("REPORTING MUTATION ANALYSIS", flush=True)
            print("-----------------", flush=True)

            mut2report.report(wk_dir, mut_prot, feature_name)

    print("\n-----------------", flush=True)
    print("FINISH", flush=True)
    print("-----------------", flush=True)


def version():
    return "1.0.1"


def run():
    global usage

    usage = "readmapper [-1 fastq_R1_.fastq] [-2 fastq_R2_.fastq] [-wd work directory] [-i " \
            "initial of the user] <-F Overwrite output directory (Default=False)> "

    parser = argparse.ArgumentParser(
        prog='mutAnalysis',
        usage=usage,
        description='mutAnalysis: count and report specific mutation - Version ' + version(),
    )

    parser.add_argument('-1', '--R1', dest="reads_1", default='', help="Reads file R1")
    parser.add_argument('-2', '--R2', dest="reads_2", default='', help="Reads file R2")
    parser.add_argument('-wd', '--wkDir', dest="workDir", default='',
                        help="Working directory")
    parser.add_argument('-i', '--initial', dest="initial", default='',
                        help="Initial of user")
    parser.add_argument('-f', '--force', dest="force", default='False',
                        help="Overwrite output directory")
    parser.add_argument('-v', '--verbose', dest="verbose", default="0",
                        help="log process to file. Options are 0 or 1  (default = 0 for no logging)")
    parser.add_argument('-V', '--version', action='version', version='parse_rep_detection-' + version(),
                        help="Prints version number")

    args = parser.parse_args()
    main(args)


if __name__ == '__main__':
    run()
