from __future__ import print_function
import os, urllib, numpy as np
import sys
import os
import desc.sltimer
timer = desc.sltimer.SLTimer()
webdir = 'https://raw.githubusercontent.com/COSMOGRAIL/PyCS/master/demo/demo1/data/'
rdbfile = 'trialcurves.txt'

url = os.path.join(webdir, rdbfile)
timer.download(url, and_read=True)
timer.computeLikelihood_MC(nsample=100,nprocess=10,rangeList=[[-10,0],[-40,-30],[-80,-60]],outName="testTutorial")
