# ==============================================================================
# License info here?
# ==============================================================================
import os

class SLTimer(object):
    '''
    Simple class for ingesting strong lens light curve data, and measuring the
    time delays.
    '''
    
    def __init__(self):
        return

    def download(self,url):
        datafile = url.split('/')[-1]
        if not os.path.isfile(datafile):
            urllib.urlretrieve(url, datafile)
        print 'Downloaded datafile:', url
        return datafile

    def read(self,lcfile):
        print("No IO enabled yet.")
        self.lc = None
        return

    def run(self,algorithm=None):
        if algorithm is None:
            print("No algorithms coded yet.")
        return
   
    def displaylightcurve(self,datafile):
        import pycs
        lcs = [
               pycs.gen.lc.rdbimport(datafile, 'A', 'mag_A', 'magerr_A', "Trial"),
               pycs.gen.lc.rdbimport(datafile, 'B', 'mag_B', 'magerr_B', "Trial"),
               pycs.gen.lc.rdbimport(datafile, 'C', 'mag_C', 'magerr_C', "Trial"),
               pycs.gen.lc.rdbimport(datafile, 'D', 'mag_D', 'magerr_D', "Trial") ]
        pycs.gen.mrg.colourise(lcs)
        lcs[1].shifttime(-5.0)
        lcs[2].shifttime(-20.0)
        lcs[3].shifttime(-70.0)
        pycs.gen.lc.display(lcs, figsize=(20, 7), jdrange=(53900, 55500))
        print "And here's your first plot!"
        pycs.gen.util.writepickle(lcs, "trialcurves_true-shifted.pkl")
        lcs = pycs.gen.util.readpickle("trialcurves_true-shifted.pkl")
        for l in lcs:
            l.resetshifts()
        return

    def splinemodel(self,datafile):
        import pycs
        lcs = [
               pycs.gen.lc.rdbimport(datafile, 'A', 'mag_A', 'magerr_A', "Trial"),
               pycs.gen.lc.rdbimport(datafile, 'B', 'mag_B', 'magerr_B', "Trial"),
               pycs.gen.lc.rdbimport(datafile, 'C', 'mag_C', 'magerr_C', "Trial"),
               pycs.gen.lc.rdbimport(datafile, 'D', 'mag_D', 'magerr_D', "Trial") ]
        pycs.gen.mrg.colourise(lcs)
        def spl(lcs):
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=50)
            for l in lcs:
                l.resetml()
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=30)
            spline = pycs.spl.topopt.opt_fine(lcs, nit=10, knotstep=20)
            return spline
        spline = spl(lcs)
        basic_time_delays = pycs.gen.lc.getnicetimedelays(lcs, separator="\n", sorted=True)
        print("Time Delays (no microlensing):")
        print(basic_time_delays)
        pycs.gen.lc.display(lcs, [spline], knotsize=0.01, figsize=(20, 7), jdrange=(53900, 55500))
        pycs.gen.polyml.addtolc(lcs[1], nparams=2, autoseasonsgap=600.0)
        pycs.gen.polyml.addtolc(lcs[2], nparams=3, autoseasonsgap=600.0)
        pycs.gen.polyml.addtolc(lcs[3], nparams=3, autoseasonsgap=600.0)
        spline = spl(lcs)
        pycs.gen.lc.display(lcs, [spline], knotsize=0.01, figsize=(20, 7), jdrange=(53900, 55500))
        polynomial_microlensing_time_delays = pycs.gen.lc.getnicetimedelays(lcs, separator="\n", sorted=True)
        print("Time Delays (microlensing included, with polynomials):")
        print(polynomial_microlensing_time_delays)
        print("cf. Time Delays (no microlensing):")
        print(basic_time_delays)
        pycs.gen.splml.addtolc(lcs[0], knotstep=150)
        pycs.gen.splml.addtolc(lcs[1], knotstep=150)
        pycs.gen.splml.addtolc(lcs[2], knotstep=150)
        pycs.gen.splml.addtolc(lcs[3], knotstep=150)
        spline = spl(lcs)
        spline_microlensing_time_delays = pycs.gen.lc.getnicetimedelays(lcs, separator="\n", sorted=True)
        print("Time Delays (microlensing included, with splines):")
        print(spline_microlensing_time_delays)
        print("cf. Time Delays (microlensing included, with polynomials):")
        print(polynomial_microlensing_time_delays)
        pycs.gen.lc.display(lcs, [spline], knotsize=0.01, figsize=(20, 7), jdrange=(53900, 55500))
        pycs.gen.util.writepickle((lcs, spline), "optspline.pkl")
        return

# ==============================================================================
