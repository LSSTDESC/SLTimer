# SLTimer
[![Build Status](https://travis-ci.org/LSSTDESC/SLTimer.svg?branch=master)](https://travis-ci.org/LSSTDESC/SLTimer)[![Coverage Status](https://coveralls.io/repos/github/LSSTDESC/SLTimer/badge.svg?branch=master)](https://coveralls.io/github/LSSTDESC/SLTimer?branch=master)

Time delay measurement from LSST light curves.

In the first instance (DC1), SLTimer will be a simple python class for calling community-developed time delay estimation algorithms, to enable early testing of DM-produced light curves in the [Twinkles project](https://github.com/LSSTDESC/Twinkles), and to check the Evil Team's simulated light curves in the [2nd Time Delay Challenge](http://timdelaychallenge.org). 

http://sltimer.readthedocs.io/

<!-- ![Example mock multi-filter LSST data from the 2nd Time Delay Challenge](https://raw.githubusercontent.com/DarkEnergyScienceCollaboration/SLTimeDelayChallenge/master/docs/TDC2/figs/example_lcs.png?token=AArY9_ujo_lqsLQGHlmkO47i7eQAMgpvks5XXu5qwA%3D%3D) -->

## Demos

Run the SLTimer IPython notebooks with [Binder](http://mybinder.org):

[![Binder](http://mybinder.org/badge.svg)](http://mybinder.org/repo/LSSTDESC/SLTimer)

Browse them on GitHub:
* [PyCS Tutorial](https://github.com/LSSTDESC/SLTimer/blob/master/notebooks/PyCS_Tutorial.ipynb): all of the steps in the `PyCS` demo1 tutorial, in notebook form.
* [SLTimer Tutorial](https://github.com/LSSTDESC/SLTimer/blob/master/notebooks/SLTimer_Tutorial.ipynb): carry out a simple "free-knot spline" `PyCS` model fit, on the `PyCS` demo1 trialcurves dataset.

## People

SLTimer development:
* [Phil Marshall](https://github.com/LSSTDESC/SLTimer/issues/new?body=@drphilmarshall) (SLAC)
* [Chun-Hao To](https://github.com/LSSTDESC/SLTimer/issues/new?body=@chto) (Stanford)
* [Milan M. Williams](https://github.com/LSSTDESC/SLTimer/issues/new?body=@milanwilliams) (Notre Dame HS, San Jose)
* [Tom Glanzman](https://github.com/LSSTDESC/SLTimer/issues/new?body=@TomGlanzman) (SLAC)

Open source time delay measurement software:
* Malte Tewes (PyCS)

## License, Contributing etc

This is open source software, available under the BSD license. If you are interested in this project, please do drop us a line via the hyperlinked contact names above, or by [writing us an issue](https://github.com/LSSTDESC/SLTimer/issues/new). To get started contributing to the SLTimer project, fork the repo and work through the [installation notes](https://github.com/LSSTDESC/SLTimer/blob/master/INSTALL.md). Pull requests are always welcome! We'd love to run your time delay measurement code :-)

