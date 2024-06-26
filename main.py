import FileMngr
from MotionTree import MotionTreeTest
from MotionTree_2 import MotionTree2
from MotionTree_Init import MotionTreeInit
from difflib import SequenceMatcher
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt


def main():
    files_dict = FileMngr.read_file_paths()
    # Read parameter file to get parameters ( Protein PDB file names, protein chains, window size, domain size, ratio )
    param_dict = FileMngr.read_param_file()
    # 10j7 = 390
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
    engine = MotionTreeTest(files_dict, param_dict)
    print(engine.protein_1.main_atoms_coords.shape)
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

