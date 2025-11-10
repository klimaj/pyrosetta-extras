"""
PyRosettaCluster environment reproducibility unit test suite using the `unittest` framework.
"""

__author__ = "Jason C. Klima"


import json
import os
import pyrosetta
import pyrosetta.distributed
import pyrosetta.distributed.io as io
import shlex
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest
import uuid

from pyrosetta.distributed.cluster import recreate_environment
from pyrosetta.distributed.cluster.config import get_environment_var


class TestEnvironmentReproducibility(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        pyrosetta.distributed.init(
            options="-run:constant_seed 1 -multithreading:total_threads 1",
            extra_options="-out:level 300 -ignore_unrecognized_res 1 -load_PDB_components 0",
            set_logging_handler="logging",
        )
        cls.workdir = tempfile.TemporaryDirectory()
        cls.run_tag = uuid.uuid4().hex[:12]

    @classmethod
    def tearDownClass(cls):
        try:
            cls.workdir.cleanup()
        except Exception as ex:
            print(f"Warning: failed to cleanup temporary directory: {ex}... Continuing.")
        os.environ.pop(get_environment_var(), None)

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

    @staticmethod
    def run_subprocess(cmd, module_dir=None, cwd=None):
        print("Running command:", cmd)
        if module_dir:
            env = os.environ.copy()
            pythonpath = os.environ.get("PYTHONPATH")
            env["PYTHONPATH"] = f"{module_dir}{os.pathsep + pythonpath if pythonpath else ''}"
        else:
            env = None
        try:
            # Use live output streaming for GitHub Actions visibility
            process = subprocess.Popen(
                shlex.split(cmd),
                cwd=cwd,
                env=env,
                shell=False,
                stdout=sys.stdout,
                stderr=sys.stderr,
                text=True,
            )
            returncode = process.wait()
            if returncode != 0:
                raise subprocess.CalledProcessError(returncode, cmd)
        except subprocess.CalledProcessError as ex:
            print(f"Subprocess command failed (return code: {ex.returncode}): {cmd}", flush=True)
            raise
        except Exception as ex:
            print(f"Unexpected error in subprocess: {ex}", flush=True)
            raise
        else:
            print(f"Return code: {returncode}", flush=True)

        return returncode

    def recreate_environment_test(self, environment_manager="conda"):
        """Test for PyRosettaCluster decoy reproducibility in a recreated virtual environment."""
        self.assertIn(environment_manager, ("conda", "mamba", "uv", "pixi"))

        test_script = os.path.join(os.path.dirname(__file__), "recreate_environment_test_runs.py")

        # Create new environment
        original_env_name = f"{environment_manager}_env_{self.run_tag}"
        original_env_dir = os.path.join(self.workdir.name, original_env_name)
        setup_env_script = os.path.join(os.path.dirname(__file__), "setup_envs.py")
        module = os.path.splitext(os.path.basename(setup_env_script))[0]
        cmd = "{0} -m {1} --env_manager '{2}' --env_dir '{3}'".format(
            sys.executable,
            module,
            environment_manager,
            original_env_dir,
        )
        returncode = TestEnvironmentReproducibility.run_subprocess(
            cmd,
            module_dir=os.path.dirname(setup_env_script),
            cwd=None,
        )
        self.assertEqual(returncode, 0, msg=f"Subprocess command failed: {cmd}")

        # Run original simulation inside new environment
        original_output_path = os.path.join(original_env_dir, f"{environment_manager}_original_outputs")
        original_scorefile_name = "test_scores.json"
        if environment_manager == "pixi":
            cmd = "pixi run python -u {0} --env_manager '{1}' --output_path '{2}' --scorefile_name '{3}'".format(
                test_script,
                environment_manager,
                original_output_path,
                original_scorefile_name,
            )
        elif environment_manager == "uv":
            cmd = "uv run -p {0} python -u {1} --env_manager '{2}' --output_path '{3}' --scorefile_name '{4}'".format(
                original_env_dir,
                test_script,
                environment_manager,
                original_output_path,
                original_scorefile_name,
            )
        elif environment_manager in ("conda", "mamba"):
            cmd = "conda run -p {0} python -u {1} --env_manager '{2}' --output_path '{3}' --scorefile_name '{4}'".format(
                original_env_dir,
                test_script,
                environment_manager,
                original_output_path,
                original_scorefile_name,
            )
        returncode = TestEnvironmentReproducibility.run_subprocess(
            cmd,
            module_dir=None,
            # For pixi, activate the original pixi environment context
            # For conda/mamba/uv, run from environment directory for consistency with pixi workflow
            cwd=original_env_dir,
        )
        self.assertEqual(returncode, 0, msg=f"Subprocess command failed: {cmd}")

        # Recreate new environment from output scorefile
        original_scorefile_path = os.path.join(original_output_path, original_scorefile_name)
        self.assertTrue(os.path.isfile(original_scorefile_path), msg=f"Missing original output scorefile: {original_scorefile_path}")
        with open(original_scorefile_path, "r") as f:
            original_data = [json.loads(line) for line in f]
        self.assertEqual(len(original_data), 1)
        original_record = original_data[0]
        self.assertIn("environment_manager", original_record["metadata"])
        self.assertEqual(original_record["metadata"]["environment_manager"], environment_manager)
        self.assertIn("decoy_name", original_record["metadata"])
        original_decoy_name = original_record["metadata"]["decoy_name"]
        # Set environment manager
        os.environ[get_environment_var()] = environment_manager
        # Recreate environment
        reproduce_env_name = f"{original_env_name}_reproduce"
        recreate_environment(
            environment_name=reproduce_env_name,
            input_file=None,
            scorefile=original_scorefile_path,
            decoy_name=original_decoy_name,
            timeout=999,
            base_dir=self.workdir.name,
        )
        reproduce_env_dir = os.path.join(self.workdir.name, reproduce_env_name)
        self.assertTrue(
            os.path.isdir(reproduce_env_dir),
            f"Reproduced '{environment_manager}' environment directory was not created: '{reproduce_env_dir}'",
        )
        if environment_manager == "uv":
            # The recreated uv environment uses the PyPI 'pyrosetta-installer' package, which does not allow specifying PyRosetta version.
            # Therefore, installing the correct PyRosetta version in the recreated uv environment depends fortuitously on a prompt
            # uv environment recreation after the original uv environment creation.
            print("Running PyRosetta installer in recreated uv environment...")
            # Run PyRosetta installer with mirror fallback
            install_script = textwrap.dedent("""
                import pyrosetta_installer
                try:
                    pyrosetta_installer.install_pyrosetta(
                        distributed=False,
                        serialization=True,
                        skip_if_installed=True,
                        mirror=0
                    )
                except Exception as e:
                    print(f"Recreated PyRosetta installation with 'mirror=0' failed: {e}. Retrying with 'mirror=1'.")
                    pyrosetta_installer.install_pyrosetta(
                        distributed=False,
                        serialization=True,
                        skip_if_installed=True,
                        mirror=1
                    )
            """)
            subprocess.run(
                ["uv", "run", "-p", str(reproduce_env_dir), "python", "-c", install_script],
                check=True,
            )

        # Run reproduction simulation inside recreated environment
        reproduce_output_path = os.path.join(reproduce_env_dir, f"{environment_manager}_reproduce_outputs")
        reproduce_scorefile_name = "test_scores.json"
        if environment_manager == "pixi":
            cmd = (
                f"pixi run python -u {test_script} "
                f"--env_manager '{environment_manager}' "
                f"--output_path '{reproduce_output_path}' "
                f"--scorefile_name '{reproduce_scorefile_name}' "
                f"--original_scorefile '{original_scorefile_path}' "
                f"--original_decoy_name '{original_decoy_name}' "
                "--reproduce"
            )
        elif environment_manager == "uv":
            cmd = (
                f"uv run -p {reproduce_env_dir} python -u {test_script} "
                f"--env_manager '{environment_manager}' "
                f"--output_path '{reproduce_output_path}' "
                f"--scorefile_name '{reproduce_scorefile_name}' "
                f"--original_scorefile '{original_scorefile_path}' "
                f"--original_decoy_name '{original_decoy_name}' "
                "--reproduce"
            )
        elif environment_manager in ("conda", "mamba"):
            cmd = (
                f"conda run -p {reproduce_env_dir} python -u {test_script} "
                f"--env_manager '{environment_manager}' "
                f"--output_path '{reproduce_output_path}' "
                f"--scorefile_name '{reproduce_scorefile_name}' "
                f"--original_scorefile '{original_scorefile_path}' "
                f"--original_decoy_name '{original_decoy_name}' "
                "--reproduce"
            )
        returncode = TestEnvironmentReproducibility.run_subprocess(
            cmd,
            module_dir=None,
            # For pixi, activate the recreated pixi environment context
            # For conda/mamba/uv, run from recreated environment directory for consistency with pixi workflow
            cwd=reproduce_env_dir,
        )
        self.assertEqual(returncode, 0, msg=f"Subprocess command failed: {cmd}")

        # Validate reproduced decoy is identical to original decoy
        reproduce_scorefile_path = os.path.join(reproduce_output_path, reproduce_scorefile_name)
        self.assertTrue(os.path.isfile(reproduce_scorefile_path), msg=f"Missing reproduced output scorefile: {reproduce_scorefile_path}")
        with open(reproduce_scorefile_path, "r") as f:
            reproduce_data = [json.loads(line) for line in f]
        self.assertEqual(len(reproduce_data), 1)
        reproduce_record = reproduce_data[0]
        self.assertIn("environment_manager", reproduce_record["metadata"])
        self.assertEqual(reproduce_record["metadata"]["environment_manager"], environment_manager)

        self.assertEqual(
            original_record["scores"]["SEQUENCE"],
            reproduce_record["scores"]["SEQUENCE"],
        )
        self.assertEqual(
            original_record["scores"]["VALUE"],
            reproduce_record["scores"]["VALUE"],
        )
        self.assertEqual(
            original_record["scores"]["total_score"],
            reproduce_record["scores"]["total_score"],
        )
        self.assertListEqual(
            original_record["instance"]["seeds"],
            reproduce_record["instance"]["seeds"],
        )
        self.assertListEqual(
            original_record["instance"]["decoy_ids"],
            reproduce_record["instance"]["decoy_ids"],
        )
        self.assertNotEqual(
            original_record["metadata"]["author"],
            reproduce_record["metadata"]["author"],
        )
        self.assertNotEqual(
            original_record["metadata"]["decoy_name"],
            reproduce_record["metadata"]["decoy_name"],
        )
        original_pose = io.pose_from_file(original_record["metadata"]["output_file"])
        reproduce_pose = io.pose_from_file(reproduce_record["metadata"]["output_file"])
        self.assert_atom_coordinates(original_pose, reproduce_pose)

    @unittest.skipIf(shutil.which("conda") is None, "The executable 'conda' is not available.")
    def test_recreate_environment_conda(self):
        return self.recreate_environment_test(environment_manager="conda")

    @unittest.skipIf(shutil.which("mamba") is None, "The executable 'mamba' is not available.")
    def test_recreate_environment_mamba(self):
        return self.recreate_environment_test(environment_manager="mamba")

    @unittest.skipIf(shutil.which("uv") is None, "The executable 'uv' is not available.")
    def test_recreate_environment_uv(self):
        return self.recreate_environment_test(environment_manager="uv")

    @unittest.skipIf(shutil.which("pixi") is None, "The executable 'pixi' is not available.")
    def test_recreate_environment_pixi(self):
        return self.recreate_environment_test(environment_manager="pixi")
