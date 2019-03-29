import unittest
import sys

import pandas as pd

from os import path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))

import pathway_assessor as _


class TestPathwayAssessor(unittest.TestCase):

    def setUp(self):
        self.expression_table_f = 'tests/expression_table_na.tsv'
        self.expression_table = _.expression_table(self.expression_table_f)

    # sanity check
    def test_hello_world_returns_str(self):
        self.assertIsInstance('hello world', str)

    def test_expression_table_is_dataframe(self):
        self.assertIsInstance(self.expression_table, pd.DataFrame)

if __name__ == '__main__':
    unittest.main()
