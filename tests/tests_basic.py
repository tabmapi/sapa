import unittest
import pytest
from unittest.mock import patch
import pathlib as pl

import sys
sys.path.append("../bin")
import sapa as s
import os
import shutil as sh

struc = s.sapa("../docs/batio3/batio3_iso.cif")
struc2 = s.sapa("../docs/scf3/scf3_alt_isodistort_v2.cif")

class TestCaseBase(unittest.TestCase):
    def assertIsFile(self, path):
        if not pl.Path(path).resolve().is_file():
            raise AssertionError(f"File does not exist: {path}")
        
    def assertIsNotFile(self, path):
        if pl.Path(path).resolve().is_file():
            raise AssertionError(f"File exists: {path}")    


class SapaTests(TestCaseBase):
    def test_find_tc_exe(self):
        for root, dirs, files in os.walk("C:\\", topdown=False):
            for name in files:
                if name == "tc.exe":
                    topaspath = os.path.join(root, name)
                    dirs[:] = []
                    break
        self.assertIn("tc.exe", topaspath)

    def test_irrep_list_batio3(self):
        irrep_list = struc.irrep_list()
        irrep_list_true = ["GM4-", "GM5-", "R5+", "R2-", "R3-", "R4-", "R5-", "X1+", "X2+", "X5+", "X3-", "X5-", "M1+", "M2+", "M3+", "M4+", "M5+", "M2-", "M5-"]
        for irrep in irrep_list_true:
            self.assertIn(irrep, irrep_list)

    def test_irrep_list_scf3(self):
        irrep_list = struc2.irrep_list()
        irrep_list_true = ["GM4-", "GM5-", "R5+", "R2-", "R3-", "R4-", "R5-", "X1+", "X2+", "X5+", "X3-", "X5-", "M1+", "M2+", "M3+", "M4+", "M5+", "M2-", "M5-"]
        for irrep in irrep_list_true:
            self.assertIn(irrep, irrep_list)
    
    def test_lprms(self):
        self.assertEqual(struc.cell_length_a, "8.01940")
        self.assertEqual(struc.cell_length_b, "8.01940")
        self.assertEqual(struc.cell_length_c, "8.01940")
        self.assertEqual(struc2.cell_length_a, "8.02490")
        self.assertEqual(struc2.cell_length_b, "8.02490")
        self.assertEqual(struc2.cell_length_c, "8.02490")


    def test_write_inp(self):
        struc2.write_inp(sample='ScF3', filenameformat = 'ScF3_##n##K_dofr.xye',
               isneutron = False,
               qmax = 21.0,
               dq = 0.08,
               startx = 1.0,
               finishx = 10.0,
               lattice = 'C',
               cycles=1)
        path = pl.Path("./batch_modes.inp")
        self.assertIsFile(path)
        os.remove(path)


    def test_execute(self):
        sh.copyfile("../docs/scf3/ScF3_125K_dofr.xye", "./ScF3_125K_dofr.xye")
        struc2.write_inp(sample='ScF3', filenameformat = 'ScF3_##n##K_dofr.xye',
               isneutron = False,
               qmax = 21.0,
               dq = 0.08,
               startx = 1.0,
               finishx = 10.0,
               lattice = 'C',
               cycles=10)
        struc2.execute(["125"])
        for irrep in struc2.irrep_list():
            path = pl.Path(f"./{irrep}_{struc2.sample}_125_out.txt")
            self.assertIsFile(path)
            with open(path, 'rb') as f:
                num_lines = sum(1 for _ in f)
                self.assertEqual(12, num_lines)
        struc2.create_hdf5(["125"])
        struc2.cleanup()
        path = pl.Path(f"./{struc2.sample}.hdf5")
        self.assertIsFile(path)
        for irrep in struc2.irrep_list():
            path = pl.Path(f"./{irrep}_{struc2.sample}_125_out.txt")
            self.assertIsNotFile(path)
        os.remove("./ScF3_125K_dofr.xye")
        os.remove("./batch_modes.out")
        os.remove(f"./{struc2.sample}.hdf5")
        os.remove(f"./{struc2.sample}.zip")
        os.remove(f"./{struc2.sample}_monitor.txt")

        
        



if __name__ == '__main__':
    unittest.main(verbosity=2)
        