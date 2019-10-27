import unittest
import sys

import numpy as np
import pandas as pd

import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pathway_assessor as _


class TestPathwayAssessor(unittest.TestCase):

    def setUp(self):

        self.test_dir = os.path.dirname(os.path.abspath(__file__))
        self.ascending = True
        self.expression_table_f = '{}/expression_table.tsv'.format(self.test_dir)
        self.expression_table_unprocessed = pd.read_csv(self.expression_table_f, sep='\t', header=0, index_col=0)

        self.expression_table = _.processed_expression_table(self.expression_table_unprocessed)

        self.expression_ranks = _.expression_ranks(self.expression_table, ascending=self.ascending)
        self.bg_genes = _.bg_genes(self.expression_ranks)

        self.sample_order = self.expression_table.columns

        self.sample_pathway = ['SLC2A6', 'PHOSPHO1', 'PIKFYVE', 'VHL']

        self.pathway_ranks = _.pathway_ranks(self.sample_pathway, self.expression_ranks, rank_method='max')

        self.effective_pathway = _.effective_pathway(self.pathway_ranks)
        self.b = _.b(self.expression_ranks, self.pathway_ranks)
        self.c = _.c(self.effective_pathway, self.pathway_ranks)
        self.d = _.d(self.bg_genes, self.pathway_ranks, self.b, self.c)

        self.sample_2x2 = _.sample_2x2(
            self.pathway_ranks.to_dict(),
            self.b.to_dict(),
            self.c.to_dict(),
            self.d.to_dict()
        )
        self.p_values = _.p_values(self.sample_2x2)

        self.harmonic_averages = self.p_values.apply(_.harmonic_average)
        self.log_harmonic_averages = _.neg_log(self.harmonic_averages)
        self.geometric_averages = self.p_values.apply(_.geometric_average)
        self.log_geometric_averages = _.neg_log(self.geometric_averages)
        self.log_min_p_vals = _.neg_log(self.p_values.min())

        self.user_pathway_f = '{}/user_pathways.txt'.format(self.test_dir)
        self.user_pathway_db, self.user_pw_data = _.user_pathways(self.user_pathway_f)

        # run with db pathway
        self.db_pathway = _.db_pathways_dict('kegg')

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

    def test_expression_ranks_is_dataframe_of_expected_values_with_max_rank_method(self):
        self.assertEqual(self.expression_ranks['Sample_A'].loc['A1'], 100)

    def test_expression_ranks_is_dataframe_of_expected_values_with_min_rank_method(self):
        expression_ranks_min_method = _.expression_ranks(self.expression_table,
                                                         ascending=self.ascending,
                                                         rank_method='min')
        self.assertEqual(expression_ranks_min_method['Sample_A'].loc['A1'], 6)

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

    def test_c_returns_df_of_expected_values(self):
        self.assertEqual(self.c['Sample_A'].loc['PIKFYVE'], 0)
        self.assertEqual(self.c['Sample_B'].loc['PIKFYVE'], 0)
        self.assertEqual(self.c['Sample_C'].loc['PIKFYVE'], 0)
        self.assertEqual(self.c['Sample_A'].loc['SLC2A6'], 1)
        self.assertEqual(self.c['Sample_B'].loc['SLC2A6'], 2)
        self.assertEqual(self.c['Sample_C'].loc['SLC2A6'], 1)
        self.assertEqual(self.c['Sample_A'].loc['PHOSPHO1'], 2)
        self.assertEqual(self.c['Sample_B'].loc['PHOSPHO1'], 1)
        self.assertTrue(pd.isna(self.c['Sample_C'].loc['PHOSPHO1']))

    def test_d_returns_df_of_expected_values(self):
        self.assertEqual(self.d['Sample_A'].loc['SLC2A6'], 95)
        self.assertEqual(self.d['Sample_A'].loc['PHOSPHO1'], 95)
        self.assertEqual(self.d['Sample_A'].loc['PIKFYVE'], 95)
        self.assertEqual(self.d['Sample_B'].loc['SLC2A6'], 96)
        self.assertEqual(self.d['Sample_B'].loc['PHOSPHO1'], 96)
        self.assertEqual(self.d['Sample_B'].loc['PIKFYVE'], 95)
        self.assertEqual(self.d['Sample_C'].loc['SLC2A6'], 95)
        self.assertTrue(pd.isna(self.d['Sample_C'].loc['PHOSPHO1']))
        self.assertEqual(self.d['Sample_C'].loc['PIKFYVE'], 95)

    def test_sample_2x2_returns_dataframe_of_numpy_arrays(self):
        self.assertIsInstance(self.sample_2x2, pd.DataFrame)
        self.assertTrue(all([type(i) == np.ndarray for i in self.sample_2x2.values]))

    def test_sample_2x2_returns_dataframe_of_expected_values(self):
        expected_dict = {
            'Sample_A': {
                'SLC2A6': [[2, 2], [1, 95]],
                'PHOSPHO1': [[1, 2], [2, 95]],
                'PIKFYVE': [[3, 2], [0, 95]],
            },
            'Sample_B': {
                'SLC2A6': [[1, 0], [2, 96]],
                'PHOSPHO1': [[2, 0], [1, 96]],
                'PIKFYVE': [[3, 1], [0, 95]],
            },
            'Sample_C': {
                'SLC2A6': [[1, 2], [1, 95]],
                'PIKFYVE': [[2, 2], [0, 95]],
            },
        }
        self.assertEqual(self.sample_2x2['Sample_A'].to_dict(), expected_dict['Sample_A'])
        self.assertEqual(self.sample_2x2['Sample_B'].to_dict(), expected_dict['Sample_B'])
        self.assertEqual(self.sample_2x2['Sample_C'].to_dict()['SLC2A6'], expected_dict['Sample_C']['SLC2A6'])
        self.assertEqual(self.sample_2x2['Sample_C'].to_dict()['PIKFYVE'], expected_dict['Sample_C']['PIKFYVE'])
        self.assertTrue(np.isnan(self.sample_2x2['Sample_C'].to_dict()['PHOSPHO1']).all())


    def test_p_values_returns_dict_of_expected_p_values(self):
        expected_dict = {
            'Sample_A': {
                'PHOSPHO1': 0.0881880024737169,
                'PIKFYVE': 6.184291898577608e-05,
                'SLC2A6': 0.0035868893011750088
            },
            'Sample_B': {
                'SLC2A6': 0.03030303030302651,
                'PHOSPHO1': 0.0006184291898577372,
                'PIKFYVE': 0.000025502234633308575,
            },
            'Sample_C': {
                'SLC2A6': 0.05998763141620062,
                'PIKFYVE': 0.0012368583797154767,
            },
        }
        sample_b = self.p_values['Sample_B'].to_dict()
        sample_c = self.p_values['Sample_C'].to_dict()
        self.assertEqual(self.p_values['Sample_A'].to_dict(), expected_dict['Sample_A'])
        self.assertAlmostEqual(sample_b['SLC2A6'], expected_dict['Sample_B']['SLC2A6'])
        self.assertAlmostEqual(sample_b['PHOSPHO1'], expected_dict['Sample_B']['PHOSPHO1'])
        self.assertAlmostEqual(sample_b['PIKFYVE'], expected_dict['Sample_B']['PIKFYVE'])
        self.assertAlmostEqual(sample_c['SLC2A6'], expected_dict['Sample_C']['SLC2A6'])
        self.assertAlmostEqual(sample_c['PIKFYVE'], expected_dict['Sample_C']['PIKFYVE'])
        self.assertTrue(pd.isna(sample_c['PHOSPHO1']))

    def test_harmonic_average_returns_expected_val(self):
        p_vals = [0.1, 0.2, 0.3, np.nan]
        expected = 0.16363636363636364
        self.assertEqual(
            _.harmonic_average(p_vals),
            expected
        )

    def test_harmonic_average_returns_zero_if_zero_p_val(self):
        p_vals = [0.1, 0.2, 0.3, 0, np.nan]
        expected = 0
        self.assertEqual(
            _.harmonic_average(p_vals),
            expected
        )

    def test_harmonic_averages_returns_dict_of_expected_values(self):
        expected_dict = {
            'Sample_A': 0.00018225855699385572,
            'Sample_B': 7.341739625203858e-05,
            'Sample_C': 0.002423742683482853
        }
        self.assertAlmostEqual(self.harmonic_averages['Sample_A'], expected_dict['Sample_A'])
        self.assertAlmostEqual(self.harmonic_averages['Sample_B'], expected_dict['Sample_B'])
        self.assertAlmostEqual(self.harmonic_averages['Sample_C'], expected_dict['Sample_C'])

    def test_log_harmonic_averages_returns_dict_of_expected_values(self):
        expected_dict = {
            'Sample_A': 8.610084236222475,
            'Sample_B': 9.519349644266763,
            'Sample_C': 6.02244237008846
        }
        self.assertAlmostEqual(self.log_harmonic_averages['Sample_A'], expected_dict['Sample_A'])
        self.assertAlmostEqual(self.log_harmonic_averages['Sample_B'], expected_dict['Sample_B'])
        self.assertAlmostEqual(self.log_harmonic_averages['Sample_C'], expected_dict['Sample_C'])

    def test_geometric_average_returns_expected_val(self):
        p_vals = [0.1, 0.2, 0.3, np.nan]
        expected = 0.181712059283214
        self.assertAlmostEqual(
            _.geometric_average(p_vals),
            expected
        )

    def test_geometric_average_returns_nan_if_all_vals_are_nan(self):
        p_vals = [np.nan, np.nan]
        self.assertTrue(
            np.isnan(_.geometric_average(p_vals))
        )

    def test_geometric_averages_returns_dict_of_expected_values(self):
        expected_dict = {
            'Sample_A': 0.002694464626458197,
            'Sample_B': 0.0007818403721029415,
            'Sample_C': 0.008613721878282994
        }
        self.assertAlmostEqual(self.geometric_averages['Sample_A'], expected_dict['Sample_A'])
        self.assertAlmostEqual(self.geometric_averages['Sample_B'], expected_dict['Sample_B'])
        self.assertAlmostEqual(self.geometric_averages['Sample_C'], expected_dict['Sample_C'])

    def test_log_geometric_averages_returns_dict_of_expected_values(self):
        expected_dict = {
            'Sample_A': 5.916555748731008,
            'Sample_B': 7.153859966001466,
            'Sample_C': 4.75439878004548
        }
        self.assertAlmostEqual(self.log_geometric_averages['Sample_A'], expected_dict['Sample_A'])
        self.assertAlmostEqual(self.log_geometric_averages['Sample_B'], expected_dict['Sample_B'])
        self.assertAlmostEqual(self.log_geometric_averages['Sample_C'], expected_dict['Sample_C'])

    def test_log_min_p_vals_returns_dict_of_expected_values(self):
        expected_dict = {
            'Sample_A': 9.690912952571226,
            'Sample_B': 10.576744476960645,
            'Sample_C': 6.695180679017199
        }
        self.assertAlmostEqual(self.log_min_p_vals['Sample_A'], expected_dict['Sample_A'])
        self.assertAlmostEqual(self.log_min_p_vals['Sample_B'], expected_dict['Sample_B'])
        self.assertAlmostEqual(self.log_min_p_vals['Sample_C'], expected_dict['Sample_C'])

    def test_user_pathway_returns_dict_of_expected_shape(self):
        expected_pw_db = {
            'EMT_kircUp': {'TGFB1', 'LAMA5', 'EGFR'},
            'EMT_kircDwn': {'CDH1', 'EGF', 'ERBB2'},
            'EMT_kirc': {'VCAN', 'EMP3', 'ZEB1', 'CDH2', 'GEM', 'COL4A2', 'SPOCK1'},
            'Sample_pathway': {'SLC2A6', 'PHOSPHO1', 'PIKFYVE', 'VHL'}
        }
        expected_pw_data = {
            'EMT_kircUp': {
                'db': 'ScientificReport',
                'count': 3,
            },
            'EMT_kircDwn': {
                'db': 'ScientificReport',
                'count': 3,
            },
            'EMT_kirc': {
                'db': 'Review',
                'count': 7
            },
            'Sample_pathway': {
                'db': 'Test_Suite',
                'count': 4
            }
        }
        self.assertIsInstance(self.user_pathway_db, dict)
        self.assertDictEqual(self.user_pathway_db, expected_pw_db)
        self.assertIsInstance(self.user_pw_data, dict)
        self.assertDictEqual(self.user_pw_data, expected_pw_data)

    def test_pathway_assessor_returns_three_dataframes_of_expected_values(self):
        user_pathways = {
            'Sample_pathway': ['SLC2A6', 'PHOSPHO1', 'PIKFYVE', 'VHL']
        }
        results = _.all(
            expression_table=self.expression_table,
            pathways=user_pathways
        )
        harmonic_avg = results['harmonic'].loc['Sample_pathway'].to_dict()
        geometric_avg = results['geometric'].loc['Sample_pathway'].to_dict()
        min_p_vals = results['min_p_val'].loc['Sample_pathway'].to_dict()
        expected_harmonic_avg = {
            'Sample_A': 8.610084236222475,
            'Sample_B': 9.519349644266763,
            'Sample_C': 6.02244237008846
        }
        expected_geometric_avg = {
            'Sample_A': 5.916555748731008,
            'Sample_B': 7.153859966001466,
            'Sample_C': 4.75439878004548
        }
        expected_min_p_vals = {
            'Sample_A': 9.690912952571226,
            'Sample_B': 10.576744476960645,
            'Sample_C': 6.695180679017199
        }
        self.assertEqual(len(results), 3)
        self.assertAlmostEqual(harmonic_avg['Sample_A'], expected_harmonic_avg['Sample_A'])
        self.assertAlmostEqual(harmonic_avg['Sample_B'], expected_harmonic_avg['Sample_B'])
        self.assertAlmostEqual(harmonic_avg['Sample_C'], expected_harmonic_avg['Sample_C'])
        self.assertAlmostEqual(geometric_avg['Sample_A'], expected_geometric_avg['Sample_A'])
        self.assertAlmostEqual(geometric_avg['Sample_B'], expected_geometric_avg['Sample_B'])
        self.assertAlmostEqual(geometric_avg['Sample_C'], expected_geometric_avg['Sample_C'])
        self.assertAlmostEqual(min_p_vals['Sample_A'], expected_min_p_vals['Sample_A'])
        self.assertAlmostEqual(min_p_vals['Sample_B'], expected_min_p_vals['Sample_B'])
        self.assertAlmostEqual(min_p_vals['Sample_C'], expected_min_p_vals['Sample_C'])

    def test_pathway_assessor_returns_none_if_statistic_is_set_to_false(self):
        results = _.all(
            expression_table=self.expression_table,
            pathways=self.user_pathway_db,
            geometric=False,
            min_p_val=False,
        )
        self.assertIsNone(results['geometric'])
        self.assertIsNone(results['min_p_val'])

    def test_db_pathways_returns_dict_of_expected_length(self):
        expected_len = 369
        self.assertIsInstance(self.db_pathway, dict)
        self.assertEqual(len(self.db_pathway), expected_len)

    def test_run_with_kegg_returns_dicts_of_dataframes_with_shape(self):
        results = _.all(
            expression_table=self.expression_table,
        )
        self.assertIsInstance(results, dict)
        self.assertIsInstance(results['harmonic'], pd.DataFrame)
        self.assertEqual(results['harmonic'].shape, (369, 3))
        self.assertEqual(results['geometric'].shape, (369, 3))
        self.assertEqual(results['min_p_val'].shape, (369, 3))

    def test_validate_pathway_returns_true_if_pathways_are_valid(self):
        valid_pw = {'pathway': ['a', 'b', 'c']}
        self.assertTrue(_.validate_pathways(valid_pw))

    def test_validate_pathway_raises_TypeError_if_pathways_are_invalid(self):
        self.assertRaises(TypeError, _.validate_pathways({'broken': 'should_be_list'}))
        self.assertRaises(TypeError, _.validate_pathways({
            'okay': ['a', 'b', 'c'],
            'broken': 'should_be_list'
        }))
        with self.assertRaises(TypeError):
            _.validate_pathways(['Should', 'be', 'dict'])

    def test_db_pathways_dict_opens_all_pathway_dicts(self):
        self.assertIsInstance(_.db_pathways_dict('hallmark'), dict)
        self.assertIsInstance(_.db_pathways_dict('reactome'), dict)
        self.assertIsInstance(_.db_pathways_dict('hmdb_smpdb'), dict)

    def test_validate_db_name_returns_true_if_pathways_are_valid(self):
        self.assertTrue(_.validate_db_name('kegg'))
        self.assertTrue(_.validate_db_name('hallmark'))
        self.assertTrue(_.validate_db_name('reactome'))
        self.assertTrue(_.validate_db_name('hmdb_smpdb'))

    def test_validate_raises_error_if_db_not_available(self):
        with self.assertRaises(ValueError):
            _.validate_db_name('nonexistent')

    def test_pa_stats_returns_expected_results_for_default_mode_harmonic(self):
        user_pathways = {
            'Sample_pathway': ['SLC2A6', 'PHOSPHO1', 'PIKFYVE', 'VHL']
        }
        results = _.pa_stats(
            expression_table=self.expression_table,
            pathways=user_pathways
        )
        results = results.loc['Sample_pathway']
        expected_results = {
            'Sample_A': 8.610084236222475,
            'Sample_B': 9.519349644266763,
            'Sample_C': 6.02244237008846
        }
        self.assertAlmostEqual(results['Sample_A'], expected_results['Sample_A'])
        self.assertAlmostEqual(results['Sample_B'], expected_results['Sample_B'])
        self.assertAlmostEqual(results['Sample_C'], expected_results['Sample_C'])

    def test_harmonic_returns_expected_results(self):
        user_pathways = {
            'Sample_pathway': ['SLC2A6', 'PHOSPHO1', 'PIKFYVE', 'VHL']
        }
        results = _.harmonic(
            expression_table=self.expression_table,
            pathways=user_pathways
        )
        results = results.loc['Sample_pathway']
        expected_results = {
            'Sample_A': 8.610084236222475,
            'Sample_B': 9.519349644266763,
            'Sample_C': 6.02244237008846
        }
        self.assertAlmostEqual(results['Sample_A'], expected_results['Sample_A'])
        self.assertAlmostEqual(results['Sample_B'], expected_results['Sample_B'])
        self.assertAlmostEqual(results['Sample_C'], expected_results['Sample_C'])

    def test_geometric_returns_expected_results(self):
        user_pathways = {
            'Sample_pathway': ['SLC2A6', 'PHOSPHO1', 'PIKFYVE', 'VHL']
        }
        results = _.geometric(
            expression_table=self.expression_table,
            pathways=user_pathways
        )
        results = results.loc['Sample_pathway']
        expected_results = {
            'Sample_A': 5.916555748731008,
            'Sample_B': 7.153859966001466,
            'Sample_C': 4.75439878004548
        }
        self.assertAlmostEqual(results['Sample_A'], expected_results['Sample_A'])
        self.assertAlmostEqual(results['Sample_B'], expected_results['Sample_B'])
        self.assertAlmostEqual(results['Sample_C'], expected_results['Sample_C'])

    def test_min_p_val_returns_expected_results(self):
        user_pathways = {
            'Sample_pathway': ['SLC2A6', 'PHOSPHO1', 'PIKFYVE', 'VHL']
        }
        results = _.min_p_val(
            expression_table=self.expression_table,
            pathways=user_pathways
        )
        results = results.loc['Sample_pathway']
        expected_results = {
            'Sample_A': 9.690912952571226,
            'Sample_B': 10.576744476960645,
            'Sample_C': 6.695180679017199
        }
        self.assertAlmostEqual(results['Sample_A'], expected_results['Sample_A'])
        self.assertAlmostEqual(results['Sample_B'], expected_results['Sample_B'])
        self.assertAlmostEqual(results['Sample_C'], expected_results['Sample_C'])


if __name__ == '__main__':
    unittest.main()
