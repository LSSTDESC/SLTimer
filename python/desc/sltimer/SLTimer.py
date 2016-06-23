# ==============================================================================
# License info here?
# ==============================================================================
import os, urllib

class SLTimer(object):
    '''
    Simple class for ingesting strong lens light curve data, and measuring the
    time delays.
    '''
    
    def __init__(self):
        return

    def download(self,url,format='pycs'):
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
    #========================Demo Model
    def display_light_curve(self):
        datafile = 'trialcurves.txt'
        import pycs
        lcs = [
               pycs.gen.lc.rdbimport(datafile, 'A', 'mag_A', 'magerr_A', "Trial"),
               pycs.gen.lc.rdbimport(datafile, 'B', 'mag_B', 'magerr_B', "Trial"),
               pycs.gen.lc.rdbimport(datafile, 'C', 'mag_C', 'magerr_C', "Trial"),
               pycs.gen.lc.rdbimport(datafile, 'D', 'mag_D', 'magerr_D', "Trial") ] #create another to accept the TDC
        pycs.gen.mrg.colourise(lcs)
        lcs[1].shifttime(-5.0)
        lcs[2].shifttime(-20.0)
        lcs[3].shifttime(-70.0)
        pycs.gen.lc.display(lcs, figsize=(20, 7), jdrange=(53900, 55500))
        print "And here's your first plot!"
        for l in lcs:
            l.resetshifts()
        pycs.gen.lc.display(lcs, filename="fig_trialcurves_true-shifted.pdf")
        return
    #==========================================Time Delays
    def basic_time_delays(self):
        import pycs
        datafile = 'trialcurves.txt'
        lcs =[
              pycs.gen.lc.rdbimport(datafile, 'A', 'mag_A', 'magerr_A', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'B', 'mag_B', 'magerr_B', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'C', 'mag_C', 'magerr_C', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'D', 'mag_D', 'magerr_D', "Trial") ]
        for l in lcs: l.resetshifts()
        for l in lcs: print(l.longinfo())
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
        return
    
    def polynomial_time_delays(self):
        import pycs
        datafile = 'trialcurves.txt'
        lcs =[
              pycs.gen.lc.rdbimport(datafile, 'A', 'mag_A', 'magerr_A', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'B', 'mag_B', 'magerr_B', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'C', 'mag_C', 'magerr_C', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'D', 'mag_D', 'magerr_D', "Trial") ]
        def spl(lcs):
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=50)
            for l in lcs:
                l.resetml()
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=30)
            spline = pycs.spl.topopt.opt_fine(lcs, nit=10, knotstep=20)
            return spline
        pycs.gen.polyml.addtolc(lcs[1], nparams=2, autoseasonsgap=600.0)
        pycs.gen.polyml.addtolc(lcs[2], nparams=3, autoseasonsgap=600.0)
        pycs.gen.polyml.addtolc(lcs[3], nparams=3, autoseasonsgap=600.0)
        spline = spl(lcs)
        polynomial_microlensing_time_delays = pycs.gen.lc.getnicetimedelays(lcs, separator="\n", sorted=True)
        print("Time Delays (microlensing included, with polynomials):")
        print(polynomial_microlensing_time_delays)
        return
    
    def spline_time_delays(self):
        import pycs
        datafile = 'trialcurves.txt'
        lcs =[
              pycs.gen.lc.rdbimport(datafile, 'A', 'mag_A', 'magerr_A', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'B', 'mag_B', 'magerr_B', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'C', 'mag_C', 'magerr_C', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'D', 'mag_D', 'magerr_D', "Trial") ]
        def spl(lcs):
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=50)
            for l in lcs:
                l.resetml()
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=30)
            spline = pycs.spl.topopt.opt_fine(lcs, nit=10, knotstep=20)
            return spline
        pycs.gen.splml.addtolc(lcs[0], knotstep=150)
        pycs.gen.splml.addtolc(lcs[1], knotstep=150)
        pycs.gen.splml.addtolc(lcs[2], knotstep=150)
        pycs.gen.splml.addtolc(lcs[3], knotstep=150)
        spline = spl(lcs)
        spline_microlensing_time_delays = pycs.gen.lc.getnicetimedelays(lcs, separator="\n", sorted=True)
        print("Time Delays (microlensing included, with splines):")
        print(spline_microlensing_time_delays)
        return

    def estimate_time_delays(self, method='pycs-spline'):
        S = SLTimer()
        if 'pycs-spline':
            S.spline_time_delays()
        elif 'pycs-polynomial':
            S.polynomial_time_delays()
        elif 'pycs-basic':
            S.basic_time_delays()
        else:
            print "Model not available"
    #=============================================Creates the basic model version
    def create_basic_model(self):
        import pycs
        S = SLTimer()
        datafile = 'trialcurves.txt'
        lcs =[
          pycs.gen.lc.rdbimport(datafile, 'A', 'mag_A', 'magerr_A', "Trial"),
          pycs.gen.lc.rdbimport(datafile, 'B', 'mag_B', 'magerr_B', "Trial"),
          pycs.gen.lc.rdbimport(datafile, 'C', 'mag_C', 'magerr_C', "Trial"),
          pycs.gen.lc.rdbimport(datafile, 'D', 'mag_D', 'magerr_D', "Trial") ]
        def spl(lcs):
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=50)
            for l in lcs:
                l.resetml()
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=30)
            spline = pycs.spl.topopt.opt_fine(lcs, nit=10, knotstep=20)
            return spline
        spline = spl(lcs)
        S.basic_time_delays()
        pycs.gen.mrg.colourise(lcs)
        pycs.gen.lc.display(lcs, [spline], knotsize=0.01, figsize=(20, 7), jdrange=(53900, 55500))
        return

       #===============================================Creates the polynomial version
    def create_polynomial_model(self):
        import pycs
        S = SLTimer()
        datafile = 'trialcurves.txt'
        lcs =[
              pycs.gen.lc.rdbimport(datafile, 'A', 'mag_A', 'magerr_A', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'B', 'mag_B', 'magerr_B', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'C', 'mag_C', 'magerr_C', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'D', 'mag_D', 'magerr_D', "Trial") ]

        def spl(lcs):
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=50)
            for l in lcs:
                l.resetml()
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=30)
            spline = pycs.spl.topopt.opt_fine(lcs, nit=10, knotstep=20)
            return spline
        pycs.gen.polyml.addtolc(lcs[1], nparams=2, autoseasonsgap=600.0)
        pycs.gen.polyml.addtolc(lcs[2], nparams=3, autoseasonsgap=600.0)
        pycs.gen.polyml.addtolc(lcs[3], nparams=3, autoseasonsgap=600.0)
        spline = spl(lcs)
        S.polynomial_time_delays()
        pycs.gen.mrg.colourise(lcs)
        pycs.gen.lc.display(lcs, [spline], knotsize=0.01, figsize=(20, 7), jdrange=(53900, 55500))
        return

    #==================================================Creates the spline model version
    def create_spline_model(self):
        import pycs
        S = SLTimer()
        datafile = 'trialcurves.txt'
        lcs =[
              pycs.gen.lc.rdbimport(datafile, 'A', 'mag_A', 'magerr_A', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'B', 'mag_B', 'magerr_B', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'C', 'mag_C', 'magerr_C', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'D', 'mag_D', 'magerr_D', "Trial") ]
        def spl(lcs):
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=50)
            for l in lcs:
                l.resetml()
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=30)
            spline = pycs.spl.topopt.opt_fine(lcs, nit=10, knotstep=20)
            return spline
        pycs.gen.splml.addtolc(lcs[0], knotstep=150)
        pycs.gen.splml.addtolc(lcs[1], knotstep=150)
        pycs.gen.splml.addtolc(lcs[2], knotstep=150)
        pycs.gen.splml.addtolc(lcs[3], knotstep=150)
        spline = spl(lcs)
        S.estimate_time_delays(method='pycs-spline')
        pycs.gen.mrg.colourise(lcs)
        pycs.gen.lc.display(lcs, [spline], knotsize=0.01, figsize=(20, 7), jdrange=(53900, 55500)) #displays spline graph
        pycs.gen.util.writepickle((lcs, spline), "optspline.pkl")
        return

    def delete_old_files(self):
        import subprocess
        subprocess.call('rm -rfv sims_copies sims_mocks', shell=True)
        subprocess.call('rm -rfv sims_copies_opt_spl sims_copies_opt_disp sims_copies_opt_regdiff', shell=True)
        subprocess.call('rm -rfv sims_mocks_opt_spl sims_mocks_opt_disp sims_mocks_opt_regdiff', shell=True)
        print "The old files have been deleted."
        return
    #=====================================================Time Delay Uncertainties

    #===================================================== Copying the Data
    def make_plain_copies(self):
        import pycs
        datafile = 'trialcurves.txt'
        lcs =[
              pycs.gen.lc.rdbimport(datafile, 'A', 'mag_A', 'magerr_A', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'B', 'mag_B', 'magerr_B', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'C', 'mag_C', 'magerr_C', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'D', 'mag_D', 'magerr_D', "Trial") ]
        n, npkl = 1, 4
        Ncopies = n*npkl
        print("Making",Ncopies,"copies of the original dataset:")
        pycs.sim.draw.multidraw(lcs, onlycopy=True, n=n, npkl=npkl, simset="copies")
        return

    def make_mock_light_curves(self):
        import pycs
        datafile = 'trialcurves.txt'
        lcs =[
              pycs.gen.lc.rdbimport(datafile, 'A', 'mag_A', 'magerr_A', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'B', 'mag_B', 'magerr_B', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'C', 'mag_C', 'magerr_C', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'D', 'mag_D', 'magerr_D', "Trial") ]
        (modellcs, modelspline)  = pycs.gen.util.readpickle("optspline.pkl")
        def Atweakml(lcs):
            return pycs.sim.twk.tweakml(lcs, beta=-1.5, sigma=0.25, fmin=1/500.0, fmax=None, psplot=False)
        def Btweakml(lcs):
            return pycs.sim.twk.tweakml(lcs, beta=-1.0, sigma=0.9, fmin=1/500.0, fmax=None, psplot=False)
        def Ctweakml(lcs):
            return pycs.sim.twk.tweakml(lcs, beta=-1.0, sigma=1.5, fmin=1/500.0, fmax=None, psplot=False)
        def Dtweakml(lcs):
            return pycs.sim.twk.tweakml(lcs, beta=-0.0, sigma=4.5, fmin=1/500.0, fmax=None, psplot=False)
        n, npkl = 1, 4
        Nmocks = n*npkl
        truetsr = 8.0
        print("Making",Nmocks,"synthetic datasets, varying time delays by +/-",truetsr/2.0,"days")
        pycs.sim.draw.saveresiduals(modellcs, modelspline)
        pycs.sim.draw.multidraw(modellcs, modelspline, n=n, npkl=npkl, simset="mocks",
                truetsr=truetsr, tweakml=[Atweakml, Btweakml, Ctweakml, Dtweakml])
        return
    #=====================================================Making Model Fits

    def make_plain_copies_model_fits(self):
        import pycs
        datafile = 'trialcurves.txt'
        lcs =[
              pycs.gen.lc.rdbimport(datafile, 'A', 'mag_A', 'magerr_A', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'B', 'mag_B', 'magerr_B', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'C', 'mag_C', 'magerr_C', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'D', 'mag_D', 'magerr_D', "Trial") ]
        def spl(lcs):
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=50)
            for l in lcs:
                l.resetml()
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=30)
            spline = pycs.spl.topopt.opt_fine(lcs, nit=10, knotstep=20)
            return spline
        pycs.sim.run.multirun("copies", lcs, spl, optset="spl", tsrand=10.0, keepopt=True)
        return

    def make_mock_light_curves_model_fits(self):
        import pycs
        datafile = 'trialcurves.txt'
        lcs =[
              pycs.gen.lc.rdbimport(datafile, 'A', 'mag_A', 'magerr_A', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'B', 'mag_B', 'magerr_B', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'C', 'mag_C', 'magerr_C', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'D', 'mag_D', 'magerr_D', "Trial") ]
        def spl(lcs):
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=50)
            for l in lcs:
                l.resetml()
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=30)
            spline = pycs.spl.topopt.opt_fine(lcs, nit=10, knotstep=20)
            return spline
        tsrand = 1.0
        pycs.sim.run.multirun("mocks", lcs, spl, optset="spl", tsrand=tsrand, keepopt=True)
        return

    def make_plain_copies_model(self): #The histogram will give the instrinic variance
        import pycs
        datafile = 'trialcurves.txt'
        lcs =[
              pycs.gen.lc.rdbimport(datafile, 'A', 'mag_A', 'magerr_A', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'B', 'mag_B', 'magerr_B', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'C', 'mag_C', 'magerr_C', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'D', 'mag_D', 'magerr_D', "Trial") ]
        def spl(lcs):
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=50)
            for l in lcs:
                l.resetml()
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=30)
            spline = pycs.spl.topopt.opt_fine(lcs, nit=10, knotstep=20)
            return spline
        dataresults = [
                pycs.sim.run.collect("sims_copies_opt_spl", "blue", "Free-knot spline technique")]
        pycs.sim.plot.hists(dataresults, r=5.0, nbins=100, showqs=False,
                filename="fig_intrinsicvariance.pdf", dataout=True)
        return
    #=====================================================Error Analysis

    def error_summary(self):
        import pycs
        datafile = 'trialcurves.txt'
        lcs =[
              pycs.gen.lc.rdbimport(datafile, 'A', 'mag_A', 'magerr_A', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'B', 'mag_B', 'magerr_B', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'C', 'mag_C', 'magerr_C', "Trial"),
              pycs.gen.lc.rdbimport(datafile, 'D', 'mag_D', 'magerr_D', "Trial") ]
        def spl(lcs):
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=50)
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=30)
            spline = pycs.spl.topopt.opt_fine(lcs, nit=10, knotstep=20)
        simresults = [
                      pycs.sim.run.collect("sims_mocks_opt_spl", "blue", "Free-knot spline technique")
        ]
        pycs.sim.plot.measvstrue(simresults, errorrange=3.5, r=5.0, nbins = 10, binclip=True, binclipr=20.0, #Creates error bars
                plotpoints=False, filename="fig_measvstrue.pdf", dataout=True)
        pycs.sim.plot.covplot(simresults, filename="fig_covplot.pdf")
        spl = (pycs.gen.util.readpickle("sims_copies_opt_spl_delays.pkl"),
               pycs.gen.util.readpickle("sims_mocks_opt_spl_errorbars.pkl"))
        pycs.sim.plot.newdelayplot([spl], rplot=6.0, displaytext=True,      #Creates a summary (of error bars and relationship bewtween measurements) plot
                filename = "fig_delays.pdf", refshifts=[{"colour":"gray", "shifts":(0, -5, -20, -70)}])
        return 

    #=====================================================Complete Error Analysis
    def complete_error_analysis(self):
        S = SLTimer()
        S.make_plain_copies()
        S.make_mock_light_curves()
        S.make_plain_copies_model_fits()
        S.make_mock_light_curves_model_fits()
        S.make_plain_copies_model()
        S.error_summary()
        return


# ==============================================================================
