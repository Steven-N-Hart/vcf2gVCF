from unittest import TestCase
from utils import bed_reader as bd


class Test(TestCase):
    def test_greater_than(self):
        value1 = [1, 100, 105]
        dict1 = dict()
        self.assertDictEqual(bd.greater_than(value1, dict1), {'1': 105})

        dict1 = {"1": 85}
        self.assertDictEqual(bd.greater_than(value1, dict1), {'1': 105})

        dict1 = {"1": 185}
        with self.assertRaises(ValueError):
            bd.greater_than(value1, dict1)
        self.assertRaisesRegex(ValueError, "Your bed file is not sorted!*")


