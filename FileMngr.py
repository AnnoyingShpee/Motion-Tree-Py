import traceback
import urllib.request
import matplotlib.pyplot as plt
import numpy as np
from itertools import groupby
from operator import itemgetter
from pathlib import Path
from scipy.cluster.hierarchy import dendrogram


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
    file_1 = temp_dict["input_path"] + "/" + temp_dict["protein1"] + ".pdb"
    file_2 = temp_dict["input_path"] + "/" + temp_dict["protein2"] + ".pdb"
    print(file_1)
    path = Path(file_1)
    print(path)
    print(path.exists())
    if not path.exists():
        try:
            urllib.request.urlretrieve(
                f"https://files.rcsb.org/download/{temp_dict['protein1']}.pdb",
                file_1
            )
        except Exception as e:
            print(e)

    path = Path(file_2)
    if not path.exists():
        try:
            urllib.request.urlretrieve(
                f"https://files.rcsb.org/download/{temp_dict['protein2']}.pdb",
                file_2
            )
        except Exception as e:
            print(e)

    return temp_dict


def get_files(input_path: str, output_path: str, pdb_1: str, chain_1: str, pdb_2: str, chain_2: str):
    temp_dict = {
        "input_path": input_path,
        "output_path": output_path,
        "protein1": pdb_1.lower(),
        "protein2": pdb_2.lower(),
        "chain1id": chain_1.upper(),
        "chain2id": chain_2.upper()
    }

    file_1 = temp_dict["input_path"] + "/" + temp_dict["protein1"] + ".pdb"
    path = Path(file_1)
    if not path.exists():
        try:
            urllib.request.urlretrieve(
                f"https://files.rcsb.org/download/{temp_dict['protein1']}.pdb",
                file_1
            )
        except Exception as e:
            traceback.print_exc()
            print(e)

    file_2 = temp_dict["input_path"] + "/" + temp_dict["protein2"] + ".pdb"
    path = Path(file_2)
    if not path.exists():
        try:
            urllib.request.urlretrieve(
                f"https://files.rcsb.org/download/{temp_dict['protein2']}.pdb",
                file_2
            )
        except Exception as e:
            traceback.print_exc()
            print(e)

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
                if param_name == "magnitude" or param_name == "dissimilarity":
                    param_val = int(tokens[1])
                elif param_name == "spatial_proximity":
                    param_val = float(tokens[1])
                temp_dict[param_name] = param_val
        fr.close()
    except Exception as e:
        print(e)
        return None
    return temp_dict


def save_results(files: dict, params: dict, data, image_type):
    proteins_folder = f"{files['protein1']}_{files['chain1id']}_{files['protein2']}_{files['chain2id']}"
    params_folder = f"sp_{params['spatial_proximity']}_dis_{params['dissimilarity']}_mag_{params['magnitude']}"
    dir_path = f"{files['output_path']}/{proteins_folder}/{params_folder}"
    path = Path(dir_path)
    if not path.is_dir():
        path.mkdir(parents=True)
    if image_type == "diff_dist_mat":
        fig_1, axis_1 = plt.subplots()
        axis_1.set_title(f"{proteins_folder} Distance Difference Matrix")
        axis_1.set_xlabel("Residue Number")
        axis_1.set_ylabel("Residue Number")
        axis_1.matshow(data)
        plt.savefig(f"{dir_path}/diff_dist_mat.png")
        # Saves the difference distance numpy array into a .npy binary file
        np.save(f"{dir_path}/diff_dist_arr.npy", data)
    elif image_type == "dendrogram":
        fig_1, axis_1 = plt.subplots()
        axis_1.set_title(f"{proteins_folder} Motion Tree")
        axis_1.set_xlabel("Residue Number")
        axis_1.set_ylabel("Magnitude (Å)")
        annotated_dendrogram(
            files, params,
            data,
            truncate_mode='lastp',
            p=12,
            leaf_rotation=90.,
            leaf_font_size=12.,
            show_contracted=True,
            annotate_above=params["magnitude"]
        )
        plt.savefig(f"{dir_path}/motion_tree.png")
    plt.close()


def annotated_dendrogram(files, params, *args, **kwargs):
    """
    Dendrogram customised for better readability.
    https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/
    :param files:
    :param params:
    :param args:
    :param kwargs:
    :return:
    """
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata


def check_for_existing_motion_tree(files, params):
    proteins_folder = f"{files['protein1']}_{files['chain1id']}_{files['protein2']}_{files['chain2id']}"
    params_folder = f"sp_{params['spatial_proximity']}_dis_{params['dissimilarity']}_mag_{params['magnitude']}"
    dir_path = f"{files['output_path']}/{proteins_folder}/{params_folder}"
    diff_dist_npy_file_path = f"{dir_path}/diff_dist_arr.npy"
    diff_dist_img_file_path = f"{dir_path}/diff_dist_mat.png"
    motion_tree_path = f"{dir_path}/motion_tree.png"

    file_path_1 = Path(diff_dist_npy_file_path)
    file_path_2 = Path(diff_dist_img_file_path)
    file_path_3 = Path(motion_tree_path)

    if file_path_1.exists() and file_path_2.exists() and file_path_3.exists():
        return diff_dist_npy_file_path, diff_dist_img_file_path, motion_tree_path
    else:
        return None, None, None


def write_bending_residues_to_file(files: dict, params: dict, nodes):
    proteins_folder = f"{files['protein1']}_{files['chain1id']}_{files['protein2']}_{files['chain2id']}"
    params_folder = f"sp_{params['spatial_proximity']}_dis_{params['dissimilarity']}_mag_{params['magnitude']}"
    file_path = f"{files['output_path']}/{proteins_folder}/{params_folder}/bend_res.info"

    try:
        fw = open(file_path, "w")
        num_nodes = len(nodes)
        fw.write(f"Protein 1 = {files['protein1']} ({files['chain1id']})\n")
        fw.write(f"Protein 2 = {files['protein2']} ({files['chain2id']})\n")
        for i in range(num_nodes - 1, -1, -1):
            fw.write(f"Effective Node = {num_nodes - i}\n")
            fw.write(f"Magnitude = {round(nodes[i]['magnitude'], 2)}\n")
            large_cluster = nodes[i]["large_cluster"]
            fw.write(f"Large Domain: {str(len(large_cluster)).ljust(3, ' ')} Residues\n")
            domain_res_str = ""
            # https://stackoverflow.com/questions/2154249/identify-groups-of-consecutive-numbers-in-a-list
            for k, g in groupby(enumerate(large_cluster), lambda ix: ix[0] - ix[1]):
                if len(domain_res_str) > 0:
                    domain_res_str = domain_res_str + " , "
                group = list(map(itemgetter(1), g))
                print(group)
                print(type(group))
                if len(group) > 1:
                    domain_res_str = domain_res_str + str(group[0]).ljust(3, ' ') + " - " + str(group[-1]).ljust(3, ' ')
                else:
                    domain_res_str = domain_res_str + str(group[0]).ljust(3, ' ')
            fw.write(f"Residues: {domain_res_str}\n")

            small_cluster = nodes[i]["small_cluster"]
            fw.write(f"Small Domain: {str(len(small_cluster)).ljust(3, ' ')} Residues\n")
            domain_res_str = ""
            for k, g in groupby(enumerate(small_cluster), lambda ix: ix[0] - ix[1]):
                if len(domain_res_str) > 0:
                    domain_res_str = domain_res_str + " , "
                group = list(map(itemgetter(1), g))
                print(group)
                print(type(group))
                if len(group) > 1:
                    domain_res_str = domain_res_str + str(group[0]).ljust(3, ' ') + " - " + str(group[-1]).ljust(3, ' ')
                else:
                    domain_res_str = domain_res_str + str(group[0]).ljust(3, ' ')
            fw.write(f"Residues: {domain_res_str}\n")
            fw.write("==============================================================\n")
        fw.close()
    except Exception as e:
        traceback.print_exc()
        print(e)


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