import traceback
import urllib.request
import matplotlib.pyplot as plt
import numpy as np
from itertools import groupby
from operator import itemgetter
from pathlib import Path
from scipy.cluster.hierarchy import dendrogram
import psycopg2
import gemmi
import pickle
from copy import deepcopy


def read_file_paths():
    temp_dict = {}
    try:
        fr = open(f"./data/input/file_mngr.txt", "r")
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


def ftp_files_to_disk(input_path: str, pdb_1: str, pdb_2: str):
    file_1 = input_path + "/" + pdb_1 + ".pdb"
    path = Path(file_1)
    if not path.exists():
        try:
            urllib.request.urlretrieve(
                f"https://files.rcsb.org/download/{pdb_1}.pdb",
                file_1
            )
        except Exception as e:
            traceback.print_exc()
            print(e)
            return 1

    file_2 = input_path + "/" + pdb_2 + ".pdb"
    path = Path(file_2)
    if not path.exists():
        try:
            urllib.request.urlretrieve(
                f"https://files.rcsb.org/download/{pdb_2}.pdb",
                file_2
            )
        except Exception as e:
            traceback.print_exc()
            print(e)
            return 1

    return 0


def read_param_file():
    temp_dict = {}
    try:
        fr = open(f"./data/input/params.txt", "r")
        lines = fr.readlines()
        for line in lines:
            if not ("#" in line):
                line = line.replace("\n", "")
                line = line.replace(" ", "")
                tokens = line.split("=")
                param_name = tokens[0]
                param_val = tokens[1]
                if param_name == "magnitude" or param_name == "clust_size" or param_name == "small_node":
                    param_val = int(tokens[1])
                elif param_name == "spatial_proximity":
                    param_val = float(tokens[1])
                temp_dict[param_name] = param_val
        fr.close()
    except Exception as e:
        print(e)
        return None
    return temp_dict


def save_results_to_disk(output_path, protein_1, chain_1, protein_2, chain_2, spat_prox, small_node, clust_size, magnitude, data, image_type):
    proteins_folder = f"{protein_1}_{chain_1}_{protein_2}_{chain_2}"
    params_folder = f"sp_{spat_prox}_node_{small_node}_clust_{clust_size}_mag_{magnitude}"
    dir_path = f"{output_path}/{proteins_folder}/{params_folder}"
    path = Path(dir_path)
    dpi = 70
    if not path.is_dir():
        path.mkdir(parents=True)
    if image_type == "diff_dist_mat":
        fig_1, axis_1 = plt.subplots()
        axis_1.set_title(f"{proteins_folder}_{params_folder} Distance Difference Matrix")
        axis_1.set_xlabel("Residue Number")
        axis_1.set_ylabel("Residue Number")
        axis_1.matshow(data)
        # plt.show()
        plt.savefig(f"{dir_path}/diff_dist_mat.png", dpi=dpi)
        # Saves the difference distance numpy array into a .npy binary file
        np.save(f"{dir_path}/diff_dist_arr.npy", data)
    elif image_type == "dendrogram":
        fig_1, axis_1 = plt.subplots()
        axis_1.set_title(f"{proteins_folder}_{params_folder} Motion Tree")
        axis_1.set_xlabel("Residue Number")
        axis_1.set_ylabel("Magnitude (Å)")
        annotated_dendrogram(
            data,
            truncate_mode='lastp',
            p=12,
            leaf_rotation=90.,
            leaf_font_size=12.,
            show_contracted=True,
            annotate_above=magnitude,
            max_d=magnitude
        )
        plt.savefig(f"{dir_path}/motion_tree.png", dpi=dpi)
    plt.close()


def annotated_dendrogram(*args, **kwargs):
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


def get_motion_tree_outputs(output_path, protein_1, chain_1, protein_2, chain_2, spat_prox, small_node, clust_size, magnitude):
    proteins_folder = f"{protein_1}_{chain_1}_{protein_2}_{chain_2}"
    params_folder = f"sp_{spat_prox}_node_{small_node}_clust_{clust_size}_mag_{magnitude}"
    dir_path = f"{output_path}/{proteins_folder}/{params_folder}"
    diff_dist_npy_file_path = f"{dir_path}/diff_dist_arr.npy"
    diff_dist_img_file_path = f"{dir_path}/diff_dist_mat.png"
    motion_tree_path = f"{dir_path}/motion_tree.png"

    file_path_1 = Path(diff_dist_npy_file_path)
    file_path_2 = Path(diff_dist_img_file_path)
    file_path_3 = Path(motion_tree_path)
    print(file_path_1)
    print(file_path_2)
    print(file_path_3)

    if file_path_1.exists() and file_path_2.exists() and file_path_3.exists():
        return diff_dist_npy_file_path, diff_dist_img_file_path, motion_tree_path
    else:
        return None, None, None


def write_to_pdb(output_path, protein_1, protein_2, spat_prox, small_node, clust_size, magnitude):
    try:
        proteins_folder = f"{protein_1.code}_{protein_1.chain_param}_{protein_2.code}_{protein_2.chain_param}"
        params_folder = f"sp_{spat_prox}_node_{small_node}_clust_{clust_size}_mag_{magnitude}"
        pdb_path = f"{output_path}/{proteins_folder}/{params_folder}/{proteins_folder}.pdb"
        fw = open(pdb_path, "w")

        protein_1_polymer = protein_1.get_polymer()
        protein_2_polymer = protein_2.get_polymer()
        util_res_1 = protein_1.utilised_res_indices
        util_res_2 = protein_2.utilised_res_indices

        ptype = protein_1_polymer.check_polymer_type()
        superimpose_result = gemmi.calculate_superposition(protein_1_polymer, protein_2_polymer,
                                                           ptype, gemmi.SupSelect.CaP)
        protein_2_polymer.transform_pos_and_adp(superimpose_result.transform)

        atom_count = 1
        subchain = protein_1.chain_param
        fw.write(f"MODEL{'1'.rjust(9, ' ')}\n")
        for i in util_res_1:
            r = protein_1_polymer[i]
            res_name = r.name.rjust(3, " ")
            res_num = str(r.seqid.num).rjust(4, " ")
            for a in r:
                atom_num = str(atom_count).rjust(5, " ")
                atom_name = a.name.ljust(4, " ")
                x = str(round(a.pos.x, 3)).rjust(8, " ")
                y = str(round(a.pos.y, 3)).rjust(8, " ")
                z = str(round(a.pos.z, 3)).rjust(8, " ")
                row = f"ATOM  {atom_num} {atom_name} {res_name} {subchain}{res_num}    {x}{y}{z}\n"
                fw.write(row)
                atom_count += 1
        fw.write("ENDMDL\n")

        atom_count = 1
        subchain = protein_2.chain_param
        fw.write(f"MODEL{'2'.rjust(9, ' ')}\n")
        for i in util_res_2:
            r = protein_2_polymer[i]
            res_name = r.name.rjust(3, " ")
            res_num = str(r.seqid.num).rjust(4, " ")
            for a in r:
                atom_num = str(atom_count).rjust(5, " ")
                atom_name = a.name.ljust(4, " ")
                x = str(round(a.pos.x, 3)).rjust(8, " ")
                y = str(round(a.pos.y, 3)).rjust(8, " ")
                z = str(round(a.pos.z, 3)).rjust(8, " ")
                row = f"ATOM  {atom_num} {atom_name} {res_name} {subchain}{res_num}    {x}{y}{z}\n"
                fw.write(row)
                atom_count += 1
        fw.write("ENDMDL\n")

        return superimpose_result

    except Exception as e:
        traceback.print_exc()
        print(e)


def write_domains_to_pml(output_path, protein_1, protein_2, spat_prox, small_node, clust_size, magnitude, nodes):
    try:
        proteins_folder = f"{protein_1.code}_{protein_1.chain_param}_{protein_2.code}_{protein_2.chain_param}"
        params_folder = f"sp_{spat_prox}_node_{small_node}_clust_{clust_size}_mag_{magnitude}"
        num_nodes = len(nodes)

        large_dom_col = "[0  ,0  ,255]"
        small_dom_col = "[255,0  ,0  ]"
        non_dom_col = "[128,128,128]"
        regions = 0

        for i in range(num_nodes-1, -1, -1):
            node_num = num_nodes - i
            pml_path = f"{output_path}/{proteins_folder}/{params_folder}/node_{node_num}.pml"
            fw = open(pml_path, "w")
            fw.write(f"load {proteins_folder}.pdb, node_{node_num}\n")

            large_domain = nodes[i]["large_domain"]
            small_domain = nodes[i]["small_domain"]

            non_domain = deepcopy(nodes[i]["large_domain"])
            non_domain.extend(small_domain)

            # Colour the large domain
            large_dom_res = protein_1.get_residue_nums(large_domain)
            groups = group_continuous_num(large_dom_res)
            first_line = True
            for group in groups:
                if first_line:
                    fw.write(f"select region{regions}, node_{node_num} and resi {group[0]}-{group[-1]}\n")
                    first_line = False
                else:
                    fw.write(f"select region{regions}, region{regions} + (node_{node_num} and resi {group[0]}-{group[-1]})\n")
            fw.write(f"set_color colour{regions} = {large_dom_col}\n")
            fw.write(f"color colour{regions}, region{regions}\n")
            fw.write("deselect\n")

            # Colour the small domain
            regions += 1
            small_dom_res = protein_1.get_residue_nums(small_domain)
            groups = group_continuous_num(small_dom_res)
            first_line = True
            for group in groups:
                if first_line:
                    fw.write(f"select region{regions}, node_{node_num} and resi {group[0]}-{group[-1]}\n")
                    first_line = False
                else:
                    fw.write(f"select region{regions}, region{regions} + (node_{node_num} and resi {group[0]}-{group[-1]})\n")
            fw.write(f"set_color colour{regions} = {small_dom_col}\n")
            fw.write(f"color colour{regions}, region{regions}\n")
            fw.write("deselect\n")

            # Colour the rest that are not domains as grey
            regions += 1
            non_dom_res = protein_1.get_residue_nums(non_domain, utilised=False)
            if len(non_dom_res) > 0:
                groups = group_continuous_num(non_dom_res)
                first_line = True
                for group in groups:
                    if first_line:
                        fw.write(f"select region{regions}, node_{node_num} and resi {group[0]}-{group[-1]}\n")
                        first_line = False
                    else:
                        fw.write(f"select region{regions}, region{regions} + (node_{node_num} and resi {group[0]}-{group[-1]})\n")
                fw.write(f"set_color colour{regions} = {non_dom_col}\n")
                fw.write(f"color colour{regions}, region{regions}\n")
                fw.write("deselect\n")

            regions += 1

    except Exception as e:
        traceback.print_exc()
        print(e)


def write_info_file(output_path, protein_1, protein_2, spat_prox, small_node, clust_size, magnitude, nodes, superimpose_result):
    proteins_folder = f"{protein_1.code}_{protein_1.chain_param}_{protein_2.code}_{protein_2.chain_param}"
    params_folder = f"sp_{spat_prox}_node_{small_node}_clust_{clust_size}_mag_{magnitude}"
    file_path = f"{output_path}/{proteins_folder}/{params_folder}/domains.info"

    try:
        fw = open(file_path, "w")
        num_nodes = len(nodes)
        fw.write(f"Protein 1 = {protein_1.code} ({protein_1.chain_param})\n")
        fw.write(f"Protein 2 = {protein_2.code} ({protein_2.chain_param})\n")
        fw.write(f"Whole Protein RMSD = {superimpose_result.rmsd}\n")
        fw.write(f"Number of Effective Nodes = {num_nodes}\n\n")
        for i in range(num_nodes - 1, -1, -1):
            fw.write("==========================================================================\n")
            fw.write(f"Effective Node {num_nodes - i}\n")
            fw.write(f"Magnitude = {round(nodes[i]['magnitude'], 2)}\n")
            fw.write("--------------------------------------------------------------------------\n")
            large_domain = nodes[i]["large_domain"]
            small_domain = nodes[i]["small_domain"]
            large_size = len(large_domain)
            small_size = len(small_domain)

            fw.write(f"{protein_1.code} ({protein_1.chain_param})\n")

            fw.write(f"Large Domain: {str(large_size).ljust(3, ' ')} Residues\n")
            large_dom_res = protein_1.get_residue_nums(large_domain)
            domain_res_str = build_info_dom_res_str(large_dom_res)
            fw.write(f"Residues: {domain_res_str}\n")

            fw.write(f"Small Domain: {str(small_size).ljust(3, ' ')} Residues\n")
            small_dom_res = protein_1.get_residue_nums(small_domain)
            domain_res_str = build_info_dom_res_str(small_dom_res)
            fw.write(f"Residues: {domain_res_str}\n\n")

            fw.write(f"{protein_2.code} ({protein_2.chain_param})\n")

            fw.write(f"Large Domain: {str(large_size).ljust(3, ' ')} Residues\n")
            large_dom_res = protein_2.get_residue_nums(large_domain)
            domain_res_str = build_info_dom_res_str(large_dom_res)
            fw.write(f"Residues: {domain_res_str}\n")

            fw.write(f"Small Domain: {str(small_size).ljust(3, ' ')} Residues\n")
            small_dom_res = protein_2.get_residue_nums(small_domain)
            domain_res_str = build_info_dom_res_str(small_dom_res)
            fw.write(f"Residues: {domain_res_str}\n\n")
        fw.close()
    except Exception as e:
        traceback.print_exc()
        print(e)


def group_continuous_num(data):
    # https://stackoverflow.com/questions/2154249/identify-groups-of-consecutive-numbers-in-a-list
    for k, g in groupby(enumerate(data), lambda ix: ix[0] - ix[1]):
        yield list(map(itemgetter(1), g))


def build_info_dom_res_str(residue_nums):
    domain_res_str = ""
    groups = group_continuous_num(residue_nums)

    for group in groups:
        if len(domain_res_str) > 0:
            domain_res_str = domain_res_str + " , "
        if len(group) > 1:
            domain_res_str = domain_res_str + str(group[0]).ljust(3, ' ') + " - " + str(group[-1]).ljust(3, ' ')
        else:
            domain_res_str = domain_res_str + str(group[0]).ljust(3, ' ')

    return domain_res_str


######################################Database##############################################


conn_str = None
try:
    pw_file = open("data/input/pw.txt", "r")
    lines = pw_file.readlines()
    for line in lines:
        cleaned_line = line.replace("\n", "")
        tokens = cleaned_line.split("=")
        if conn_str is None:
            conn_str = f"{tokens[0]}='{tokens[1]}'"
        else:
            if tokens[0] == "password":
                conn_str += f" {tokens[0]}={tokens[1]}"
            else:
                conn_str += f" {tokens[0]}='{tokens[1]}'"
except Exception as e:
    print(e)

conn = None
try:
    conn = psycopg2.connect(conn_str)
    conn.autocommit = True
    cur = conn.cursor()
    cur.execute("SET SEARCH_PATH TO motion_tree_py, public;")
except Exception as e:
    print(e)


# def check_protein_exists(protein, chain):
#     try:
#         cur.execute(
#             """
#             SELECT EXISTS(SELECT 1 FROM proteins
#             WHERE protein_id=%s AND protein_chain=%s);
#             """,
#             (protein, chain)
#         )
#         is_exists = cur.fetchone()
#         return is_exists[0]
#     except Exception as e:
#         traceback.print_exc()
#         print(e)
#         return -1
#
#
# def insert_protein_dist_mat(protein, chain, dist_mat):
#     try:
#         dist_mat_bytes = pickle.dumps(dist_mat)
#         cur.execute(
#             """
#             INSERT INTO proteins (protein_id, protein_chain, dist_mat)
#             VALUES (%s, %s, %s);
#             """,
#             (protein, chain, dist_mat_bytes)
#         )
#         return 0
#     except Exception as e:
#         traceback.print_exc()
#         print(e)
#         return -1
#
#
# def get_protein_dist_mat(protein, chain):
#     # conn = psycopg2.connect(conn_str)
#     try:
#         # conn.autocommit = True
#         # cur = conn.cursor()
#         # cur.execute("SET SEARCH_PATH TO motion_tree_py, public;")
#         cur.execute(
#             """
#             SELECT protein_id, protein_chain, dist_mat FROM proteins
#             WHERE protein_id=%s AND protein_chain=%s;
#             """,
#             (protein, chain)
#         )
#         row = cur.fetchone()
#         dist_mat = pickle.loads(row[2])
#         return dist_mat
#     except Exception as e:
#         traceback.print_exc()
#         print(e)
#         return -1


def check_motion_tree_exists(protein_1, chain_1, protein_2, chain_2, spat_prox):
    try:
        cur.execute(
            """
            SELECT EXISTS(
            SELECT 1 FROM motion_trees WHERE protein_1=%s AND chain_1=%s AND protein_2=%s AND chain_2=%s AND 
            spatial_proximity=%s
            )
            """,
            (protein_1, chain_1, protein_2, chain_2, spat_prox)
        )
        row = cur.fetchone()
        return row[0]
    except Exception as e:
        traceback.print_exc()
        print(e)
        return -1


def get_motion_tree(protein_1, chain_1, protein_2, chain_2, spat_prox):
    try:
        cur.execute(
            """
            SELECT diff_dist_mat, link_mat FROM motion_trees
            WHERE protein_1=%s AND chain_1=%s AND protein_2=%s AND chain_2=%s AND spatial_proximity=%s;
            """,
            (protein_1, chain_1, protein_2, chain_2, spat_prox)
        )
        rows = cur.fetchall()
        diff_dist_mat = pickle.loads(rows[0][0])
        link_mat = pickle.loads(rows[0][1])
        return diff_dist_mat, link_mat
    except Exception as e:
        traceback.print_exc()
        print(e)
        return -1, -1


def insert_motion_tree(protein_1, chain_1, protein_2, chain_2, spat_prox, seq_identity, time_taken, diff_dist_mat, link_mat, motion_tree_exist):
    try:
        diff_dist_bin = pickle.dumps(diff_dist_mat)
        link_bin = pickle.dumps(link_mat)
        if motion_tree_exist:
            cur.execute(
                """
                UPDATE motion_trees 
                SET time_taken=%s, diff_dist_mat=%s, link_mat=%s
                WHERE protein_1=%s AND chain_1=%s AND protein_2=%s AND chain_2=%s AND spatial_proximity=%s;
                """,
                (round(time_taken, 2), diff_dist_bin, link_bin, protein_1, chain_1, protein_2, chain_2, spat_prox)
            )
        else:
            cur.execute(
                """
                INSERT INTO motion_trees
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s);
                """,
                (protein_1, chain_1, protein_2, chain_2, spat_prox, round(seq_identity, 2), round(time_taken, 2), diff_dist_bin, link_bin)
            )
        return 0
    except Exception as e:
        traceback.print_exc()
        print(e)
        return -1


def check_nodes_exist(protein_1, chain_1, protein_2, chain_2, spat_prox, small_node, clust_size, magnitude):
    try:
        cur.execute(
            """
            SELECT EXISTS(
            SELECT * FROM nodes WHERE protein_1=%s AND chain_1=%s AND protein_2=%s AND chain_2=%s AND spatial_proximity=%s
            AND small_node_size=%s AND cluster_size=%s AND magnitude=%s 
            );
            """,
            (protein_1, chain_1, protein_2, chain_2, spat_prox, small_node, clust_size, magnitude)
        )
        row = cur.fetchone()
        return row[0]
    except Exception as e:
        traceback.print_exc()
        print(e)
        return -1


def get_nodes(protein_1, chain_1, protein_2, chain_2, spat_prox, small_node, clust_size, magnitude):
    try:
        cur.execute(
            """
            SELECT node, distance, large_domain, small_domain FROM nodes
            WHERE protein_1=%s AND chain_1=%s AND protein_2=%s AND chain_2=%s AND spatial_proximity=%s AND
            small_node_size=%s AND cluster_size=%s AND magnitude=%s;
            """,
            (protein_1, chain_1, protein_2, chain_2, spat_prox, small_node, clust_size, magnitude)
        )
        # List of tuples
        rows = cur.fetchall()
        nodes = {}
        for row in rows:
            nodes[row[0]] = {
                "magnitude": row[1],
                "large_domain": [],
                "small_domain": []
            }
            large_domain = pickle.loads(row[2])
            small_domain = pickle.loads(row[3])
            for i in range(large_domain.shape[0]):
                temp = []
                for j in range(large_domain[i][0], large_domain[i][1]+1):
                    temp.append(j)
                nodes[row[0]]["large_domain"].extend(temp)
            for i in range(small_domain.shape[0]):
                temp = []
                for j in range(small_domain[i][0], small_domain[i][1]+1):
                    temp.append(j)
                nodes[row[0]]["small_domain"].extend(temp)

        return nodes
    except Exception as e:
        traceback.print_exc()
        print(e)
        return -1


def insert_nodes(protein_1, chain_1, protein_2, chain_2, spat_prox, small_node, clust_size, magnitude, nodes, nodes_exist):
    try:
        if nodes_exist:
            cur.execute(
                """
                DELETE FROM nodes
                WHERE protein_1=%s AND chain_1=%s AND protein_2=%s AND chain_2=%s AND spatial_proximity=%s 
                AND small_node_size=%s AND cluster_size=%s AND magnitude=%s;
                """,
                (protein_1, chain_1, protein_2, chain_2, spat_prox, small_node, clust_size, magnitude)
            )
        for key in nodes.keys():
            large_domain_groups = group_continuous_num(nodes[key]["large_domain"])
            temp = []
            for group in large_domain_groups:
                temp.append([group[0], group[-1]])
            large_domain_bin = pickle.dumps(np.asarray(temp))

            small_domain_groups = group_continuous_num(nodes[key]["small_domain"])
            temp = []
            for group in small_domain_groups:
                temp.append([group[0], group[-1]])
            small_domain_bin = pickle.dumps(np.asarray(temp))

            cur.execute(
                """
                INSERT INTO nodes
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);
                """,
                (protein_1, chain_1, protein_2, chain_2, spat_prox, small_node, clust_size, magnitude, key, nodes[key]["magnitude"], large_domain_bin, small_domain_bin)
            )
        return 0
    except Exception as e:
        traceback.print_exc()
        print(e)
        return -1

