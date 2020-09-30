from csv import DictReader


def read_mutation_database(mutation_file):
    mut_dict = {}
    with open(mutation_file, "r") as mut_file:
        reader = DictReader(mut_file, delimiter="\t")
        for row in reader:
            if row["Proteic"] == None:
                row["Proteic"] = ""
            if row["Nucleic"] == None:
                row["Nucleic"] = ""
            mut_dict[row["Sequence"]] = {
                "proteic": row["Proteic"].replace("[", "").replace("]", "").split(","),
                "nucleic": row["Nucleic"].replace("[", "").replace("]", "").split(",")
            }
    return mut_dict
