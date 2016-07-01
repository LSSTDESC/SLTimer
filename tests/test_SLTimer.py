"""
    Test code for SLTimer class.
    """
from __future__ import absolute_import, division
import unittest

import desc.sltimer

class SLTimerTestCase(unittest.TestCase):
    "TestCase class for SLTimer class."
    def setUp(self):
        """
            Set up each test with a new SLTimer object.
            """
        self.timer = desc.sltimer.SLTimer()

    def test_download_and_read_in(self):
        '''
            Checking to see if the URL downloaded properly:
            '''
        url = "https://raw.githubusercontent.com/COSMOGRAIL/PyCS/master/demo/demo1/data/trialcurves.txt"
        self.timer.download(url)
        numLines = 0
        with open(self.timer.datafile, 'r') as file:
            for line in file:
                wordsList = line.split()
                numLines += 1
        self.assertEqual(numLines, 194)
        return


if __name__ == '__main__':
    unittest.main()
