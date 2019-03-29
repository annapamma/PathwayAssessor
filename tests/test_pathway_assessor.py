import unittest
import sys

import pandas as pd

from os import path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))

import pathway_assessor as _


class TestPathwayAssessor(unittest.TestCase):

    def setUp(self):
        self.ascending = True
        self.expression_table_f = 'tests/expression_table.tsv'
        self.expression_table = _.expression_table(self.expression_table_f)
        self.expression_ranks = _.expression_ranks(self.expression_table, ascending=self.ascending)

    # sanity check
    def test_hello_world_returns_str(self):
        self.assertIsInstance('hello world', str)

    def test_expression_table_is_dataframe_of_expected_values(self):
        self.assertIsInstance(self.expression_table, pd.DataFrame)
        self.assertAlmostEqual(self.expression_table['Sample_A'].loc['PIKFYVE'], 7.442)
        self.assertAlmostEqual(self.expression_table['Sample_B'].loc['PIKFYVE'], 10.969)
        self.assertAlmostEqual(self.expression_table['Sample_C'].loc['PIKFYVE'], 5.520)
        self.assertTrue(pd.isna(self.expression_table['Sample_C'].loc['PHOSPHO1']))

    def test_expression_ranks_is_dataframe_of_expected_values(self):
        self.assertIsInstance(self.expression_ranks, pd.DataFrame)
        self.assertEqual(self.expression_ranks['Sample_A'].loc['PIKFYVE'], 5)
        self.assertEqual(self.expression_ranks['Sample_B'].loc['PIKFYVE'], 4)
        self.assertEqual(self.expression_ranks['Sample_C'].loc['PIKFYVE'], 4)



if __name__ == '__main__':
    unittest.main()
