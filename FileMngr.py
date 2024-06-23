import sys
from pathlib import Path
import urllib.request


def read_file_paths():
    temp_dict = {}
    try:
        fr = open(f"data/input/file_mngr.txt", "r")
        lines = fr.readlines()
        for line in lines:
            if not ("#" in line):
                line = line.replace("\n", "")
                line = line.replace(" ", "")
                tokens = line.split("=")
                param_name = tokens[0]
                param_val = tokens[1]
                temp_dict[param_name] = param_val
        fr.close()
    except Exception as e:
        print(e)

    temp_dict["protein1"] = temp_dict["protein1"].lower()
    temp_dict["protein2"] = temp_dict["protein2"].lower()
    file_1 = temp_dict["input_path"] + temp_dict["protein1"] + ".pdb"
    file_2 = temp_dict["input_path"] + temp_dict["protein2"] + ".pdb"

    path = Path(file_1)
    if not path.exists():
        try:
            urllib.request.urlretrieve(
                f"http://files.rcsb.org/download/{temp_dict['protein1']}.pdb",
                file_1
            )
        except Exception as e:
            print(e)
            sys.exit()

    path = Path(file_2)
    if not path.exists():
        try:
            urllib.request.urlretrieve(
                f"http://files.rcsb.org/download/{temp_dict['protein2']}.pdb",
                file_2
            )
        except Exception as e:
            print(e)
            sys.exit()

    return temp_dict


def read_param_file():
    temp_dict = {}
    try:
        fr = open(f"data/input/params.txt", "r")
        lines = fr.readlines()
        for line in lines:
            if not ("#" in line):
                line = line.replace("\n", "")
                line = line.replace(" ", "")
                tokens = line.split("=")
                param_name = tokens[0]
                param_val = tokens[1]
                if param_name == "threshold" or param_name == "avg_diff":
                    param_val = int(tokens[1])
                temp_dict[param_name] = param_val
        fr.close()
    except Exception as e:
        print(e)
        return None
    return temp_dict


def write_clustering(file_name: str, data, n=0):
    file_path = f"data/output/clustering/{file_name}.txt"

    try:
        if n:
            fo = open(file_path, "a")
        else:
            fo = open(file_path, "w")
        fo.write(data)
        fo.write("\n")
        fo.close()
    except Exception as e:
        print(e)
        print(data)


# def write_pymol_file(pml_file_name: str, pdb_file_name: str, data):
#     try:
#         fw = open(f"{output_pymol_file_path}{pml_file_name}.pml", "w")
#         fw.write("reinitialize\n")
#         fw.write(f"load {pdb_file_name}\n")
#         fw.write("bg_color white")
#         fw.write("color grey")
#         fw.close()
#     except Exception as e:
#         print(e)
#         return False
#     return True