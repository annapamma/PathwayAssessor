import unittest
import sys


class TestPathwayAssessor(unittest.TestCase):
    # sanity check
    def test_hello_world_returns_str(self):
        self.assertIsInstance('hello world', str)


if __name__ == '__main__':
    unittest.main()
