# Getting started with SLTimer

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
