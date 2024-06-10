import FileMngr
from MotionTree import MotionTree


def main():
    files_dict = FileMngr.read_file_paths()
    # Read parameter file to get parameters ( Protein PDB file names, protein chains, window size, domain size, ratio )
    param_dict = FileMngr.read_param_file()
    # Initialise MotionTree object
    engine = MotionTree(files_dict, param_dict)
    # Run the Engine
    engine.run()


if __name__ == '__main__':
    main()

