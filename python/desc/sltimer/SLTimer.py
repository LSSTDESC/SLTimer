# ==============================================================================
# License info here?
# ==============================================================================
import os, urllib

from desc.sltimer import *

class SLTimer(object):
    '''
    Simple class for ingesting strong lens light curve data, and measuring the
    time .
    '''
    
    def __init__(self):
        return

    def download(self,url,format='pycs'):
        self.datafile = url.split('/')[-1]
        if not os.path.isfile(self.datafile):
            urllib.urlretrieve(url, self.datafile)
        print 'Downloaded datafile:', url
        return
            
    def read_in(self,format=None):
        if format == 'rdb':  self.lcs = read_in_rdb_data(self.datafile)
        if format == 'tdc2': self.lcs = read_in_tdc2_data(self.datafile)
        return
    
    def spl(self):
        global spl
        def spl(lcs):
            import pycs
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=50)
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=30)
            spline = pycs.spl.topopt.opt_fine(lcs, nit=10, knotstep=20)
            return spline

    def basic_spl(self):
        global spl
        def spl(lcs):
            import pycs
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=50)
            for l in lcs:
                l.resetml()
            spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=30)
            spline = pycs.spl.topopt.opt_fine(lcs, nit=10, knotstep=20)
            return spline

    def defines_variables(self):
        self.defines_datafile_type()
        self.spl()
        return

    #==========================================================Demo Model
    def display_light_curves(self):
        import pycs
        self.defines_datafile_type()
        pycs.gen.mrg.colourise(lcs)
        lcs[1].shifttime(-5.0)
        lcs[2].shifttime(-20.0)
        lcs[3].shifttime(-70.0)
        pycs.gen.lc.display(lcs, figsize=(20, 7), jdrange=(53900, 55500))
        lcs = pycs.gen.util
        for l in lcs:
            l.resetshifts()
        pycs.gen.lc.display(lcs, filename="fig_trialcurves_true-shifted.pdf")
        return
    #=====================================================Time Delays

    def show_basic_time_delays(self):
        import pycs
        self.lcs()
        self.basic_spl()
        pycs.gen.mrg.colourise(lcs)
        spline = spl(lcs)
        basic_time_delays = pycs.gen.lc.getnicetimedelays(lcs, separator="\n", sorted=True)
        print("Time Delays (no microlensing):")
        print(basic_time_delays)
        return
    
    def show_polynomial_time_delays(self):
        import pycs
        self.defines_variables()
        pycs.gen.polyml.addtolc(lcs[1], nparams=2, autoseasonsgap=600.0)
        pycs.gen.polyml.addtolc(lcs[2], nparams=3, autoseasonsgap=600.0)
        pycs.gen.polyml.addtolc(lcs[3], nparams=3, autoseasonsgap=600.0)
        spline = spl(lcs)
        polynomial_microlensing_time_delays = pycs.gen.lc.getnicetimedelays(lcs, separator="\n", sorted=True)
        print("Time Delays (microlensing included, with polynomials):")
        print(polynomial_microlensing_time_delays)
        return
    
    def show_spline_time_delays(self):
        import pycs
        self.defines_variables()
        pycs.gen.splml.addtolc(lcs[0], knotstep=150)
        pycs.gen.splml.addtolc(lcs[1], knotstep=150)
        pycs.gen.splml.addtolc(lcs[2], knotstep=150)
        pycs.gen.splml.addtolc(lcs[3], knotstep=150)
        spline = spl(lcs)
        spline_microlensing_time_delays = pycs.gen.lc.getnicetimedelays(lcs, separator="\n", sorted=True)
        print("Time Delays (microlensing included, with splines):")
        print(spline_microlensing_time_delays)
        return


    def estimate_time_delays(self,method='pycs',microlensing='polynomial',agn='spline',error=None): #Provides both polynomial and spline time delays
        if 'pycs':
            print ("You are using the pycs method.")
        else:
            return
        if 'polynomial':
            self.polynomial_time_delays()
    
        else:
            return
        if 'spline':
            self.spline_time_delays()
        else:
            return
        if 'intrinsic variance':
            self.find_intrinsic_variance
        elif 'complete':
            self.estimate_uncertainties()
        else:
            return

    #=============================================Creates the basic model version
    def display_basic_method_plot(self):
        import pycs
        self.defines_variables()
        spline = spl(lcs)
        self.show_basic_time_delays()
        pycs.gen.mrg.colourise(lcs)
        pycs.gen.lc.display(lcs, [spline], knotsize=0.01, figsize=(20, 7), jdrange=(53900, 55500))
        return

       #===============================================Creates the polynomial version
    def display_polynomial_method_plot(self):
        import pycs
        self.display_basic_method_plot()
        pycs.gen.polyml.addtolc(lcs[1], nparams=2, autoseasonsgap=600.0)
        pycs.gen.polyml.addtolc(lcs[2], nparams=3, autoseasonsgap=600.0)
        pycs.gen.polyml.addtolc(lcs[3], nparams=3, autoseasonsgap=600.0)
        spline = spl(lcs)
        self.show_polynomial_time_delays()
        pycs.gen.mrg.colourise(lcs)
        pycs.gen.lc.display(lcs, [spline], knotsize=0.01, figsize=(20, 7), jdrange=(53900, 55500))
        return

    #==================================================Creates the spline model version
    def display_spline_method_plot(self):
        import pycs
        self.display_basic_method_plot()
        pycs.gen.splml.addtolc(lcs[0], knotstep=150)
        pycs.gen.splml.addtolc(lcs[1], knotstep=150)
        pycs.gen.splml.addtolc(lcs[2], knotstep=150)
        pycs.gen.splml.addtolc(lcs[3], knotstep=150)
        spline = spl(lcs)
        self.show_spline_time_delays()
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
    def make_plain_copies(self,n=None,npkl=None):
        import pycs
        self.defines_variables()
        Ncopies = n*npkl
        print("Making",Ncopies,"copies of the original dataset:")
        pycs.sim.draw.multidraw(lcs, onlycopy=True, n=n, npkl=npkl, simset="copies")
        return
    
    def make_mock_light_curves(self,n=None,npkl=None):
        import pycs
        self.defines_variables()
        (modellcs, modelspline)  = pycs.gen.util.readpickle("optspline.pkl")
        def Atweakml(lcs):
            return pycs.sim.twk.tweakml(lcs, beta=-1.5, sigma=0.25, fmin=1/500.0, fmax=None, psplot=False)
        def Btweakml(lcs):
            return pycs.sim.twk.tweakml(lcs, beta=-1.0, sigma=0.9, fmin=1/500.0, fmax=None, psplot=False)
        def Ctweakml(lcs):
            return pycs.sim.twk.tweakml(lcs, beta=-1.0, sigma=1.5, fmin=1/500.0, fmax=None, psplot=False)
        def Dtweakml(lcs):
            return pycs.sim.twk.tweakml(lcs, beta=-0.0, sigma=4.5, fmin=1/500.0, fmax=None, psplot=False)
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
        self.defines_variables()
        pycs.sim.run.multirun("copies", lcs, spl, optset="spl", tsrand=10.0, keepopt=True)
        return

    def make_mock_light_curves_model_fits(self):
        import pycs
        self.defines_variables()
        tsrand = 1.0
        pycs.sim.run.multirun("mocks", lcs, spl, optset="spl", tsrand=tsrand, keepopt=True)
        return

    def make_plain_copies_model(self): #The histogram will give the instrinic variance
        import pycs
        self.defines_variables()
        dataresults = [
                pycs.sim.run.collect("sims_copies_opt_spl", "blue", "Free-knot spline technique")]
        pycs.sim.plot.hists(dataresults, r=5.0, nbins=100, showqs=False,
                filename="fig_intrinsicvariance.pdf", dataout=True)
        return
    #=========================================================Error Analysis
    def error_summary(self):
        import pycs
        self.defines_variables()
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
    def estimate_uncertainties(self,n=None,npkl=None):
        self.make_plain_copies(n=n,npkl=npkl)
        self.make_mock_light_curves(n=n,npkl=npkl)
        self.make_plain_copies_model_fits()
        self.make_mock_light_curves_model_fits()
        self.make_plain_copies_model()
        self.error_summary()
        return
    
    def find_intrinsic_variance(self,n=None,npkl=None):
        self.make_plain_copies(n=n,npkl=npkl)
        self.make_plain_copies_model_fits()
        self.make_plain_copies_model()
        return


# ==============================================================================
