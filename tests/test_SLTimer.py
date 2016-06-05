"""
Test code for SLTimer class.
"""
from __future__ import absolute_import, division
import unittest
from desc.sltimer import SLTimer


class SLTimerTestCase(unittest.TestCase):
    "TestCase class for SLTimer class."
    def setUp(self):
        """
        Set up each test with a new SLTimer object.
        """
        self.timer = SLTimer()

    def test_read(self):
        "Test the read function."
        self.timer.read("lightcurve.txt")
        self.assertEqual(self.timer.lc, None)


if __name__ == '__main__':
    unittest.main()
