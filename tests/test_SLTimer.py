"""
Test code for SLTimer class.
"""
from __future__ import absolute_import, division
import unittest

import desc.sltimer
from desc.sltimer import SLTimer

class SLTimerTestCase(unittest.TestCase):
    "TestCase class for SLTimer class."
    def setUp(self):
        """
        Set up each test with a new SLTimer object.
        """
        self.timer = SLTimer()
    
    def test_download(self):
        '''
        Checking to see if the URL downloaded properly:
        '''
        url = "https://raw.githubusercontent.com/COSMOGRAIL/PyCS/master/demo/demo1/data/trialcurves.txt"
        datafile = self.timer.download(url)
        numLines = 0
        with open(datafile, 'r') as file:
            for line in file:
                wordsList = line.split()
                numLines += 1
        self.assertEqual(numLines, 194)
        return

    def test_read(self):
        "Testing the read function:"
        self.timer.read("lightcurve.txt")
        self.assertEqual(self.timer.lc, None)



if __name__ == '__main__':
    unittest.main()
