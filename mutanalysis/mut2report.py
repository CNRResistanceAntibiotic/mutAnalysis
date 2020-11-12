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
import csv
import os
import re
from csv import DictReader

from Bio.Seq import Seq


def report(work_dir, mutation, gene_name):
    threshold_view_mut = 10

    final_result = os.path.join(work_dir, "final_result.tsv")

    combinaison_final_dict = {}

    for file in os.listdir(work_dir):
        if "_count.csv" in file:
            path = os.path.join(work_dir, file)

            with open(path, "r") as file_r:
                reader = DictReader(file_r, delimiter="\t")
                combinaison_ref = ""
                combinaison_dict = {}
                depth = count = 0
                for row in reader:
                    combinaison_tmp_dict = {}
                    combinaison_ref += row["reference"]
                    depth += int(row["total_depth"])
                    if int(row["A_depth"]) != 0:
                        if count == 0:
                            combinaison_dict["A"] = int(row["A_depth"])
                        else:
                            for comb, depth_l in combinaison_dict.items():
                                if len(comb) == count:
                                    comb = comb + "A"
                                    if int(row["A_depth"]) >= depth_l:
                                        pass
                                    else:
                                        depth_l = int(row["A_depth"])
                                    combinaison_tmp_dict[comb] = depth_l
                    if int(row["C_depth"]) != 0:
                        if count == 0:
                            combinaison_dict["C"] = int(row["C_depth"])
                        else:
                            for comb, depth_l in combinaison_dict.items():
                                if len(comb) == count:
                                    comb = comb + "C"
                                    if int(row["C_depth"]) >= depth_l:
                                        pass
                                    else:
                                        depth_l = int(row["C_depth"])
                                    combinaison_tmp_dict[comb] = depth_l
                    if int(row["G_depth"]) != 0:
                        if count == 0:
                            combinaison_dict["G"] = int(row["G_depth"])
                        else:
                            for comb, depth_l in combinaison_dict.items():
                                if len(comb) == count:
                                    comb = comb + "G"
                                    if int(row["G_depth"]) >= depth_l:
                                        pass
                                    else:
                                        depth_l = int(row["G_depth"])
                                    combinaison_tmp_dict[comb] = depth_l
                    if int(row["T_depth"]) != 0:
                        if count == 0:
                            combinaison_dict["T"] = int(row["T_depth"])
                        else:
                            for comb, depth_l in combinaison_dict.items():
                                if len(comb) == count:
                                    comb = comb + "T"
                                    if int(row["T_depth"]) >= depth_l:
                                        pass
                                    else:
                                        depth_l = int(row["T_depth"])
                                    combinaison_tmp_dict[comb] = depth_l
                    if count != 0:
                        combinaison_dict = combinaison_tmp_dict
                    count += 1
                combinaison_final_dict = {}
                median_ref_depth = int(depth/3)
                for comb, depth in combinaison_dict.items():
                    if (depth/median_ref_depth)*100 >= threshold_view_mut:
                        combinaison_final_dict[comb] = depth
                    else:
                        continue
    if combinaison_final_dict:
        with open(final_result, "w") as output:
            writer = csv.writer(output, delimiter="\t")
            writer.writerow(["Gene", "Mutation", "Mean Depth", "Result", "Sensible/Resistant"])
            for comb, depth in combinaison_final_dict.items():
                pattern = re.compile('([a-zA-Z_-]+)*([0-9]*)([a-zA-Z_-]+)')
                match = pattern.match(mutation)
                if match:
                    acide_s = match.groups()[0]
                    acide_r = match.groups()[2]
                dna = Seq(comb)
                if dna.translate(table=11) == acide_s:
                    resu_type = "Sensible"
                elif dna.translate(table=11) == acide_r:
                    resu_type = "Resistant"
                else:
                    resu_type = "X"
                writer.writerow([gene_name, mutation, median_ref_depth, "{0}(depth:{1};ratio{2}%)".format(comb, depth, int((depth/median_ref_depth)*100)), resu_type])
