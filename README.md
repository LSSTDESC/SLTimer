# SLTimer
[![Build Status](https://travis-ci.org/DarkEnergyScienceCollaboration/SLTimer.svg?branch=master)](https://travis-ci.org/DarkEnergyScienceCollaboration/SLTimer)[![Coverage Status](https://coveralls.io/repos/github/DarkEnergyScienceCollaboration/SLTimer/badge.svg?branch=master)](https://coveralls.io/github/DarkEnergyScienceCollaboration/SLTimer?branch=master)

Time delay measurement from LSST light curves.

In the first instance (DC1), SLTimer will be a simple python class for calling community-developed time delay estimation algorithms, to enable early testing of DM-produced light curves in the [Twinkles project](https://github.com/DarkEnergyScienceCollaboration/Twinkles).

## Set-up

You'll need to install all the required packages. These include the open source ["python curve-shifting" package, `PyCS`](http://pycs.readthedocs.io/en/latest/), which is the first community-developed tool we are investigating. `PyCS` requires `pymc`, which is best installed via
```
$ conda install -c https://conda.binstar.org/pymc pymc
```
See the [`pymc` docs]() for more assistance. Once you have `pymc` installed, you can get the rest of the packages you need by doing
```
$ pip install -r requirements.txt
```
At this point it might be a good idea to [take a look at the `PyCS` tutorial](http://pycs.readthedocs.io/en/latest/tutorial/tutorial.html).

Finally, we need to get `SLTimer` onto our `PYTHONPATH`, using the `setup` script. From bash, do:
```
$ source <SLTimer install directory>/setup/setup.sh
```
From tcsh, do:
```
$ source <SLTimer install directory>/setup/setup.csh
```

## Demo

Placeholder:
```
    python sltimer/SLTimer.py
```

## People

* [Tom Glanzman](https://github.com/DarkEnergyScienceCollaboration/SLTimer/issues/new?body=@TomGlanzman) (SLAC)
* [Phil Marshall](https://github.com/DarkEnergyScienceCollaboration/SLTimer/issues/new?body=@drphilmarshall) (SLAC)

## License etc

This is open source software, available under the BSD license. If you are interested in this project, please do drop us a line via the hyperlinked contact names above, or by [writing us an issue](https://github.com/DarkEnergyScienceCollaboration/SLTimer/issues/new).
