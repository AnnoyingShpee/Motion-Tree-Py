import FileMngr
from Engine import Engine


def main():
    # Read command file to get parameters ( Protein PDB file names, protein chains, window size, domain size, ratio )
    file_paths_dict = FileMngr.read_file_paths()
    param_dict = FileMngr.read_param_file()
    # Concatenate PDB file names with path
    pdb_path_1 = file_paths_dict["file_path"] + file_paths_dict["filename1"]
    pdb_path_2 = file_paths_dict["file_path"] + file_paths_dict["filename2"]
    # Initialise Engine object
    engine = Engine(pdb_path_1, pdb_path_2, file_paths_dict, param_dict)
    # Run the Engine
    engine.run()


if __name__ == '__main__':
    main()

