"""Verify header-axis slicing without downloading the public grid."""
import importlib.util
from pathlib import Path
import re
import unittest

spec = importlib.util.spec_from_file_location('prepare', Path(__file__).with_name('prepare_apogee.py'))
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


class HeaderTests(unittest.TestCase):
    def header(self):
        block = """ &SYNTH
 N_OF_DIM=7
 N_P=5 11 6 8 15 10 11
 LLIMITS=-0.52288 -1.5 -0.5 -0.75 -2.5 0 3500
 STEPS=0.30103 0.25 0.5 0.25 0.25 0.5 250
 NPCA=627 627 617
 NPIX=900
"""
        block += '\n'.join(f" LABEL({i})='{name}'" for i,name in
                           enumerate(['micro','C','N','alpha','metal','logg','teff'],1))
        return (block+'\n /\n')*4 + 'PCA numerical data stays untouched\n'

    def test_five_dimensions(self):
        result = prepare.slice_header(self.header(),5)
        self.assertEqual(result.count('N_OF_DIM = 5'),4)
        self.assertEqual(result.count('N_P = 6 8 15 10 11'),4)
        self.assertEqual(result.count("LABEL(1) ='N'"),4)
        self.assertEqual(result.count("LABEL(5) ='teff'"),4)
        self.assertNotIn('LABEL(6)',result)
        self.assertTrue(result.endswith('PCA numerical data stays untouched\n'))
        self.assertEqual(result.count('NPCA=627 627 617'),4)

    def test_four_dimensions_compatible(self):
        result = prepare.slice_header(self.header())
        self.assertEqual(result.count('N_P = 8 15 10 11'),4)
        self.assertEqual(result.count("LABEL(1) ='alpha'"),4)

    def test_reject_unexpected_source(self):
        with self.assertRaises(ValueError):
            prepare.slice_header(self.header().replace('5 11 6 8 15 10 11','5 10 6 8 15 10 11'),5)
        with self.assertRaises(ValueError):
            prepare.slice_header(self.header(),6)


if __name__ == '__main__':
    unittest.main()
