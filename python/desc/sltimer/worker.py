# ======================================================================
# License info here?
# ======================================================================
from __future__ import absolute_import
import os, urllib
import subprocess
import pycs
from .reading import *

__all__ = ['SLTimer', 'spl']

class SLTimer(object):
    '''
    Worker class for ingesting strongly lensed image light curves, and
    measuring the time delays between them.
    '''

    def __init__(self):
        self.agn = None
        self.microlensing = None
        self.time_delays = None
        self.datafile = None
        self.lcs = None
        return

    def download(self, url, format='rdb', and_read=False):
        '''
        Downloads the datafile from a url.

        Parameters
        ----------
        url : string
            Web address of datafile.
        format : string
            Data format, 'rdb' or 'tdc2'
        and_read : boolean
            Read in data after downloading file?

        Notes
        -----
        Don't forget to set `and_read=True` if you want to use the data!
        '''
        self.datafile = url.split('/')[-1]
        if not os.path.isfile(self.datafile):
            urllib.urlretrieve(url, self.datafile)
        print 'Downloaded datafile:', url
        if and_read:
            self.read_in(format=format)
        return

    def read_in(self, datafile='self', format=None):
        '''
        Reads in light curve data from a file.
        '''
        if datafile == 'self':
            pass
        else:
            self.datafile = datafile

        if format == 'rdb':
            self.lcs = read_in_rdb_data(self.datafile)
        elif format == 'tdc2':
            self.lcs = read_in_tdc2_data(self.datafile)
        else:
            raise ValueError('Unrecognized or null format '+str(format))

        self.Nim = len(self.lcs)

        return

    def optimize_spline_model(self):
        '''
        Optimizes a spline model for the intrinsic variability.
        '''
        return spl(self.lcs)

    #========================================================== Plotting light curves

    def display_light_curves(self, filename=None, jdrange=(None)):
        '''
        Displays the lightcurves in a single panel plot.
        '''
        pycs.gen.mrg.colourise(self.lcs)
        # Replace the following with an optional input list of shifts
        #lcs[1].shifttime(-5.0)
        #lcs[2].shifttime(-20.0)
        #lcs[3].shifttime(-70.0)
        pycs.gen.lc.display(self.lcs, [self.agn], figsize=(20, 7), jdrange=jdrange)
        # lcs = pycs.gen.util
        # for l in lcs:
        #     l.resetshifts()
        if filename is not None:
            pycs.gen.lc.display(self.lcs, filename=filename)
        return

    def whiten(self):
        '''
        Whitens a set of multi-filter light curves to a single fictitious band.
        '''
        self.lcs = whiten(self.lcs)
        return

    #===================================================== Microlensing

    def add_polynomial_microlensing(self):
        '''
        Adds polynomial microlensing to each lightcurve.
        '''
        pycs.gen.polyml.addtolc(self.lcs[0], nparams=3,
                                autoseasonsgap=600.0)
        pycs.gen.polyml.addtolc(self.lcs[1], nparams=3,
                                autoseasonsgap=600.0)
        if self.Nim == 4:
            pycs.gen.polyml.addtolc(self.lcs[2], nparams=3,
                                    autoseasonsgap=600.0)
            pycs.gen.polyml.addtolc(self.lcs[3], nparams=3,
                                    autoseasonsgap=600.0)
        return

    def add_spline_microlensing(self):
        '''
        Adds spline microlensing to each light curve.
        '''
        pycs.gen.splml.addtolc(self.lcs[0], knotstep=150)
        pycs.gen.splml.addtolc(self.lcs[1], knotstep=150)
        if self.Nim == 4:
            pycs.gen.splml.addtolc(self.lcs[2], knotstep=150)
            pycs.gen.splml.addtolc(self.lcs[3], knotstep=150)
        return

    #========================================= Primary workhorse method

    def estimate_time_delays(self, method='pycs', microlensing='spline', agn='spline', error=None, quietly=False):
        '''
        Measures time delays between images, by modeling all the light
        curves.

        Parameters
        ----------
        method : string
            Modeling package to use.
        microlensing : string
            Choice of microlensing model to use.
        agn : string
            Choice of intrinsic AGN variability model to use.
        error : string
            Error estimation options [None, 'complete', 'intrinsic variance']

        Notes
        -----
        Provides both polynomial and spline time delays.

        Parameters
        ----------
        method: string
            Modeling package to use (currently only `pycs` is available)
        microlensing: string
            Model choice for microlensing light curves
        agn: string
            Model choice for intrinsic AGN variability
        error: boolean
            Estimate errors?
        quietly: boolean
            Redirect output to /dev/null?
        '''

        if method == 'pycs':
            # print "You are using the pycs method."
            pass
        else:
            print "The only available method is 'pycs' - exiting."
            return

        if quietly:
            as_requested = {'stdout':None, 'stderr':None}
        else:
            as_requested = {'stdout':sys.stdout, 'stderr':sys.stderr}

        # Tell the lightcurves that their model is going to include microlensing:
        if microlensing == 'polynomial':
            with SilentOperation(**as_requested):
                self.add_polynomial_microlensing()
        elif microlensing == 'spline':
            with SilentOperation(**as_requested):
                self.add_spline_microlensing()
        else:
            pass
        # Keep a record:
        self.microlensing = microlensing

        # Optimize the model for both microlensing and intrinsic variability:
        if agn == 'spline':
            with SilentOperation(**as_requested):
                self.agn = self.optimize_spline_model()
        else:
            print "Error: only free-knot spline models are available  for AGN variability at present."
            return

        # Do error analysis, if required:
        if error == 'complete':
            with SilentOperation(**as_requested):
                self.estimate_uncertainties()
        elif error == 'intrinsic variance':
            with SilentOperation(**as_requested):
                self.find_intrinsic_variance()
        else:
            return


    def initialize_time_delays(self, method=None, pars=None):
        '''
        Initializes the curve shifts by specifying 1 or 3 time delays.
        '''
        if method is None:
            dt = {'AB':0.0}
            if self.Nim == 4:
                dt['AC'] = 0.0
                dt['AD'] = 0.0

        elif method == 'guess':
            dt = pars
            assert pars is not None
            assert len(dt) == (self.Nim - 1)
            assert type(dt) == dict

        else:
            raise ValueError("Unrecognized initialization method '"+method+"'")

        # Set the shifts of each light curve object in lcs:
        # All lenses:
        self.lcs[1].shifttime(dt['AB'])
        # Quads only:
        if self.Nim == 4:
            self.lcs[2].shifttime(dt['AC'])
            self.lcs[3].shifttime(dt['AD'])

        # Report that shifting has occurred, and report time delays:
        print "Initialization completed, using method '"+method+"'"
        self.report_time_delays()

        return

    #===================================================== Resimulating the Data

    def delete_old_files(self):
        '''
        Deletes the old files from previous error simulations.
        '''
        subprocess.call('rm -rfv sims_copies sims_mocks', shell=True)
        subprocess.call('rm -rfv sims_copies_opt_spl sims_copies_opt_disp sims_copies_opt_regdiff', shell=True)
        subprocess.call('rm -rfv sims_mocks_opt_spl sims_mocks_opt_disp sims_mocks_opt_regdiff', shell=True)
        print "The old files have been deleted."
        return

    def make_plain_copies(self, n=None, npkl=None):
        '''
        Makes copies of the data.
        '''
        Ncopies = n*npkl
        print "Making", Ncopies, "copies of the original dataset:"
        pycs.sim.draw.multidraw(self.lcs, onlycopy=True, n=n, npkl=npkl, simset="copies")
        return

    def make_mock_light_curves(self, n=None, npkl=None):
        '''
        Make mock lightcurves to help estimate uncertainties.
        '''
        modellcs, modelspline = self.lcs, self.agn
        def Atweakml(xlcs):
            return pycs.sim.twk.tweakml(xlcs, beta=-1.5, sigma=0.25, fmin=1/500.0, fmax=None, psplot=False)
        def Btweakml(xlcs):
            return pycs.sim.twk.tweakml(xlcs, beta=-1.0, sigma=0.9, fmin=1/500.0, fmax=None, psplot=False)
        def Ctweakml(xlcs):
            return pycs.sim.twk.tweakml(xlcs, beta=-1.0, sigma=1.5, fmin=1/500.0, fmax=None, psplot=False)
        def Dtweakml(xlcs):
            return pycs.sim.twk.tweakml(xlcs, beta=-0.0, sigma=4.5, fmin=1/500.0, fmax=None, psplot=False)
        Nmocks = n*npkl
        truetsr = 8.0
        print "Making", Nmocks, "synthetic datasets, varying time delays by +/-", truetsr/2.0, "days"
        pycs.sim.draw.saveresiduals(modellcs, modelspline)
        pycs.sim.draw.multidraw(modellcs, modelspline, n=n, npkl=npkl, simset="mocks", truetsr=truetsr, tweakml=[Atweakml, Btweakml, Ctweakml, Dtweakml])
        return

    #========================================Making Multiple Model Fits

    def make_spline_model_fits_of_plain_copies(self):
        # Pass the optimizer function to multirun:
        pycs.sim.run.multirun("copies", self.lcs, spl, optset="spl", tsrand=10.0, keepopt=True)
        return

    def make_spline_model_fits_of_mock_light_curves(self):
        tsrand = 1.0
        # Pass the optimizer function to multirun:
        pycs.sim.run.multirun("mocks", self.lcs, spl, optset="spl", tsrand=tsrand, keepopt=True)
        return

    def plot_intrinsic_variance_histograms(self): #The histogram will give the instrinsic variance
        dataresults = [pycs.sim.run.collect("sims_copies_opt_spl", "blue", "Free-knot spline technique")]
        pycs.sim.plot.hists(dataresults, r=5.0, nbins=100, showqs=False,
                            filename="fig_intrinsicvariance.pdf", dataout=True)
        return

    #=================================================== Error Analysis

    def error_summary(self):
        simresults = [
                      pycs.sim.run.collect("sims_mocks_opt_spl", "blue", "Free-knot spline technique")]
        # Nice to replace self.time_delays with a version including error bars here...
        # Maybe write out the "samples" for post-processing! Could also make a corner plot...
        # Compare measured time delays with truth:
        pycs.sim.plot.measvstrue(simresults, errorrange=3.5, r=5.0, nbins = 1, binclip=True, binclipr=20.0,
                                 plotpoints=False, filename="fig_measvstrue.pdf", dataout=True)
        # Plot covariances between delays:
        pycs.sim.plot.covplot(simresults, filename="fig_covplot.pdf")
        # Create a summary plot (of error bars and relationship bewtween measurements):
        spl = (pycs.gen.util.readpickle("sims_copies_opt_spl_delays.pkl"),
               pycs.gen.util.readpickle("sims_mocks_opt_spl_errorbars.pkl"))
        # One last plot:
        pycs.sim.plot.newdelayplot([spl], rplot=6.0, displaytext=True,
                filename = "fig_delays.pdf", refshifts=[{"colour":"gray", "shifts":(0, -5, -20, -70)}])
        return

    #=====================================================Complete Error Analysis

    def estimate_uncertainties(self, n=None, npkl=None):
        self.delete_old_files()
        self.make_plain_copies(n=n, npkl=npkl)
        self.make_mock_light_curves(n=n, npkl=npkl)
        # Add in an option to use regdiff and disp here
        self.make_spline_model_fits_of_plain_copies()
        self.make_spline_model_fits_of_mock_light_curves()
        self.plot_intrinsic_variance_histograms()
        self.error_summary()
        return

    def find_intrinsic_variance(self,n=None, npkl=None):
        self.make_plain_copies(n=n, npkl=npkl)
        self.make_spline_model_fits_of_plain_copies()
        self.plot_intrinsic_variance_histograms()
        return

    def report_time_delays(self):
        print "Time Delays:"
        self.time_delays = pycs.gen.lc.getnicetimedelays(self.lcs, separator="\n", sorted=True)
        print self.time_delays
        return

# ======================================================================
# End of the SLTimer class.
# ======================================================================

# Optimizer functions (could go in "optimize.py" instead?)

def spl(lcs):
    spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=50)
    spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=30)
    spline = pycs.spl.topopt.opt_fine(lcs, nit=10, knotstep=20)
    return spline
