import sys
import DataMngr
from MotionTree import MotionTree
from MotionTree_2 import MotionTree2
from MotionTree_Init import MotionTreeInit
from difflib import SequenceMatcher
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt


def main():
    files_dict = DataMngr.read_file_paths()
    # Read parameter file to get parameters ( Protein PDB file names, protein chains, window size, domain size, ratio )
    param_dict = DataMngr.read_param_file()
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
    engine = MotionTree(files_dict, param_dict)
    engine.init_protein(1)
    engine.init_protein(2)
    engine.check_sequence_identity()
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

