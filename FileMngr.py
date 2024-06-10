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


# def write_rotation_vec_to_pdb(file_name: str, slide_window_residues, slide_window_indices, rotation_vectors):
#     """
#     Writes the rotation vectors of each residue of the slide window into a input file
#     :param file_name:
#     :param slide_window_residues:
#     :param slide_window_indices:
#     :param rotation_vectors:
#     :return:
#     """
#     try:
#         fw = open(f"{output_pdb_file_path}{file_name}_rot_vecs.input", "w")
#         start_index = slide_window_indices[0]
#         for r in range(len(slide_window_residues)):
#             residue_name = slide_window_residues[r].name
#             residue_num = str(start_index + r + 1).rjust(3, " ")
#             x = str(round(rotation_vectors[r][0], 3)).rjust(8, " ")
#             y = str(round(rotation_vectors[r][1], 3)).rjust(8, " ")
#             z = str(round(rotation_vectors[r][2], 3)).rjust(8, " ")
#             row = f"ATOM         CA  {residue_name} A {residue_num}    {x}{y}{z}\n"
#             fw.write(row)
#         fw.close()
#     except Exception as e:
#         print(e)
#         return False
#     return True
#
#
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