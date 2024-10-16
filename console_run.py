import FileMngr
from MotionTree import MotionTree


def main():
    files_dict = FileMngr.read_file_paths()
    # Read parameter file to get parameters ( Protein PDB file names, protein chains, window size, domain size, ratio )
    param_dict = FileMngr.read_param_file()
    # Initialise MotionTree object
    # print("Motion Tree Init")
    # engine = MotionTreeInit(files_dict, param_dict)
    # cluster_1 = engine.run()
    # print("===============================================================================")
    # print("Motion Tree 2")
    # engine = MotionTree2(files_dict, param_dict)
    # cluster_2 = engine.run()
    # print("===============================================================================")
    print("Motion Tree Test")
    engine = MotionTree(files_dict["input_path"], files_dict["output_path"],
                        files_dict["protein1"], files_dict["chain1id"], files_dict["protein2"], files_dict["chain2id"],
                        param_dict["spatial_proximity"], param_dict["small_node"],
                        param_dict["clust_size"], param_dict["magnitude"])
    engine.init_protein(1)
    engine.init_protein(2)
    engine.check_rmsd_standard()
    engine.check_sequence_identity_standard()
    engine.preprocessing_standard()
    engine.dist_mat_processing()
    engine.create_distance_difference_matrix()
    engine.run()
    # cluster_3, link_mat = engine.run()

    # for i in range(len(cluster_1[0])):
    #     print(i, cluster_1[0][i], cluster_2[0][i])

    # print(SequenceMatcher(a=cluster_1[0], b=cluster_2[0]).ratio())
    # print(SequenceMatcher(a=cluster_1[0], b=cluster_3[778]).ratio())

    # sch.dendrogram(link_mat, 30)
    # plt.show()


if __name__ == '__main__':
    main()

