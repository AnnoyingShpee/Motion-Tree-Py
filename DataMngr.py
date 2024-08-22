import traceback
import numpy as np
import psycopg2
import pickle
from itertools import groupby
from operator import itemgetter


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


def check_protein_pair_exists(protein_1, chain_1, protein_2, chain_2):
    try:
        cur.execute(
            """
            SELECT EXISTS(
            SELECT 1 FROM proteins WHERE protein_1=%s AND chain_1=%s AND protein_2=%s AND chain_2=%s
            )
            """,
            (protein_1, chain_1, protein_2, chain_2)
        )
        row = cur.fetchone()
        return row[0]
    except Exception as e:
        traceback.print_exc()
        print(e)
        return -1


def insert_protein_pair(protein_1, chain_1, protein_2, chain_2, rmsd, diff_dist_mat):
    try:
        diff_dist_bin = pickle.dumps(diff_dist_mat)
        cur.execute(
            """
            INSERT INTO proteins (protein_1, chain_1, protein_2, chain_2, rmsd, diff_dist_mat) 
            VALUES (%s, %s, %s, %s, %s, %s);
            """,
            (protein_1, chain_1, protein_2, chain_2, rmsd, diff_dist_bin)
        )
        return 0
    except Exception as e:
        traceback.print_exc()
        print(e)
        return -1


def get_protein_pair(protein_1, chain_1, protein_2, chain_2):
    try:
        cur.execute(
            """
            SELECT rmsd, diff_dist_mat FROM proteins 
            WHERE protein_1=%s AND chain_1=%s AND protein_2=%s AND chain_2=%s;
            """,
            (protein_1, chain_1, protein_2, chain_2)
        )
        row = cur.fetchone()
        diff_dist_mat = pickle.loads(row[1])
        return row[0], diff_dist_mat
    except Exception as e:
        traceback.print_exc()
        print(e)
        return -1, -1


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
            SELECT link_mat FROM motion_trees
            WHERE protein_1=%s AND chain_1=%s AND protein_2=%s AND chain_2=%s AND spatial_proximity=%s;
            """,
            (protein_1, chain_1, protein_2, chain_2, spat_prox)
        )
        row = cur.fetchone()
        link_mat = pickle.loads(row[0])
        return link_mat
    except Exception as e:
        traceback.print_exc()
        print(e)
        return -1


def insert_motion_tree(protein_1, chain_1, protein_2, chain_2, spat_prox, seq_identity, time_taken, link_mat, motion_tree_exist):
    try:
        link_bin = pickle.dumps(link_mat)
        if motion_tree_exist:
            cur.execute(
                """
                UPDATE motion_trees 
                SET time_taken=%s, link_mat=%s
                WHERE protein_1=%s AND chain_1=%s AND protein_2=%s AND chain_2=%s AND spatial_proximity=%s;
                """,
                (round(time_taken, 2), link_bin, protein_1, chain_1, protein_2, chain_2, spat_prox)
            )
        else:
            cur.execute(
                """
                INSERT INTO motion_trees
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s);
                """,
                (protein_1, chain_1, protein_2, chain_2, spat_prox, round(seq_identity, 2), round(time_taken, 2), link_bin)
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


def group_continuous_num(data):
    # https://stackoverflow.com/questions/2154249/identify-groups-of-consecutive-numbers-in-a-list
    for k, g in groupby(enumerate(data), lambda ix: ix[0] - ix[1]):
        yield list(map(itemgetter(1), g))

