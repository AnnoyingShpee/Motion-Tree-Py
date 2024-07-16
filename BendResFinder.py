import gemmi


class BendResFinder:
    def __init__(self, motion_tree):
        self.motion_tree = motion_tree

        # Bending residues indices
        self.bending_residues = {}

    def get_rotation_vectors(self):
        return
