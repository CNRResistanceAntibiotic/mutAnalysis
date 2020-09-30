import csv
import os
import re
from csv import DictReader

from Bio.Seq import Seq


def report(work_dir, mutation, gene_name):
    threshold_view_mut = 10

    final_result = os.path.join(work_dir, "final_result.tsv")

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

    with open(final_result, "w") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(["Gene", "Mutation", "Mean Depth", "Result", "S/R"])
        for comb, depth in combinaison_final_dict.items():
            pattern = re.compile('([a-zA-Z_-]+)*([0-9]*)([a-zA-Z_-]+)')
            match = pattern.match(mutation)
            if match:
                acide_s = match.groups()[0]
                acide_r = match.groups()[2]
            dna = Seq(comb)
            if dna.translate(table=11) == acide_s:
                resu_type = "S"
            elif dna.translate(table=11) == acide_r:
                resu_type = "R"
            else:
                resu_type = "X"
            writer.writerow([gene_name, mutation, median_ref_depth, "{0}(depth:{1};ratio{2}%)".format(comb, depth, int((depth/median_ref_depth)*100)), resu_type])
