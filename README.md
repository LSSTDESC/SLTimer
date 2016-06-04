# SLTimer

Time delay measurement from LSST light curves.

In the first instance (DC1), SLTimer will be a simple python class for calling community-developed time delay estimation algorithms, to enable early testing of DM-produced light curves in the [Twinkles project](https://github.com/DarkEnergyScienceCollaboration/Twinkles).

## Set-up

You'll need to install all the required packages:
```
$ pip install -r requirements.txt
```
These include the open source ["python curve-shifting" package, `PyCS`](http://pycs.readthedocs.io/en/latest/), which is the first community-developed tool we are investigating.

Then, from bash:
```
$ source <SLTimer install directory>/setup/setup.sh
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
