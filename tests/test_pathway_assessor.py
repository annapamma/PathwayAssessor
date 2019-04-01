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
        self.bg_genes = _.bg_genes(self.expression_ranks)

        self.sample_order = self.expression_table.columns

        self.sample_pathway = ['SLC2A6', 'PHOSPHO1', 'PIKFYVE', 'VHL']

        self.pathway_ranks = _.pathway_ranks(self.sample_pathway, self.expression_ranks)
        self.effective_pathway = _.effective_pathway(self.pathway_ranks)
        self.b = _.b(self.expression_ranks, self.pathway_ranks)


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

    def test_bg_returns_series_of_expected_lengths(self):
        expected = {
            'Sample_A': 100,
            'Sample_B': 99,
            'Sample_C': 99
        }
        self.assertIsInstance(self.bg_genes, pd.Series)
        self.assertDictEqual(self.bg_genes.to_dict(), expected)

    def test_pathway_rank_returns_expected_ranks(self):
        self.assertEqual(self.pathway_ranks['Sample_A'].loc['PIKFYVE'], 3)
        self.assertEqual(self.pathway_ranks['Sample_B'].loc['PIKFYVE'], 3)
        self.assertEqual(self.pathway_ranks['Sample_C'].loc['PIKFYVE'], 2)

    def test_effective_pathway_returns_series_of_expected_lengths(self):
        expected = {
            'Sample_A': 3,
            'Sample_B': 3,
            'Sample_C': 2
        }
        self.assertDictEqual(self.effective_pathway.to_dict(), expected)

    def test_b_returns_df_of_expected_values(self):
        self.assertEqual(self.b['Sample_A'].loc['PIKFYVE'], 2)
        self.assertEqual(self.b['Sample_B'].loc['PIKFYVE'], 1)
        self.assertEqual(self.b['Sample_C'].loc['PIKFYVE'], 2)
        self.assertEqual(self.b['Sample_A'].loc['SLC2A6'], 2)
        self.assertEqual(self.b['Sample_B'].loc['SLC2A6'], 0)
        self.assertEqual(self.b['Sample_C'].loc['SLC2A6'], 2)
        self.assertEqual(self.b['Sample_A'].loc['PHOSPHO1'], 2)
        self.assertEqual(self.b['Sample_B'].loc['PHOSPHO1'], 0)
        self.assertTrue(pd.isna(self.b['Sample_C'].loc['PHOSPHO1']))




if __name__ == '__main__':
    unittest.main()
