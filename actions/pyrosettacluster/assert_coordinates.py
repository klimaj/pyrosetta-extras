__author__ = "Jason C. Klima"


import argparse
import pyrosetta
import pyrosetta.distributed.io as io
import unittest


class TestAtomCoordinates(unittest.TestCase):
    def assert_atom_coordinates(self, pose1, pose2):
        self.assertEqual(pose1.size(), pose2.size())
        for res in range(1, pose1.size() + 1):
            res1 = pose1.residue(res)
            res2 = pose2.residue(res)
            self.assertEqual(res1.name(), res2.name())
            self.assertEqual(res1.natoms(), res2.natoms())
            for atom in range(1, res1.natoms() + 1):
                self.assertEqual(res1.atom_name(atom), res2.atom_name(atom))
                for axis in "xyz":
                    self.assertEqual(
                        float(getattr(res1.atom(atom).xyz(), axis)),
                        float(getattr(res2.atom(atom).xyz(), axis)),
                    )

    def assert_rmsd(self, pose1, pose2):
        rmsd = pyrosetta.rosetta.core.scoring.all_atom_rmsd_incl_hydrogens(pose1, pose2)
        self.assertEqual(rmsd, 0.0)

    def assert_energy(self, pose1, pose2):
        scorefxn = io.create_score_function()
        self.assertEqual(scorefxn(pose1), scorefxn(pose2))

    def test_coordinates(self):
        original_pose = io.pose_from_file(self.original_output_file).pose
        reproduce_pose = io.pose_from_file(self.reproduce_output_file).pose
        self.assert_atom_coordinates(original_pose, reproduce_pose)
        self.assert_rmsd(original_pose, reproduce_pose)
        self.assert_energy(original_pose, reproduce_pose)


if __name__ == "__main__":
    print("Running: {0}".format(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument('--original_output_file', type=str, required=True)
    parser.add_argument('--reproduce_output_file', type=str, required=True)
    args, remaining_argv = parser.parse_known_args()
    # Inject args into the class before running test
    TestAtomCoordinates.original_output_file = args.original_output_file
    TestAtomCoordinates.reproduce_output_file = args.reproduce_output_file
    # Run test
    unittest.main(argv=[__file__] + remaining_argv)
