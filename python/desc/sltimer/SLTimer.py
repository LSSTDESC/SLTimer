# ======================================================================
# License info here?
# ======================================================================
from __future__ import absolute_import
import os, urllib
import subprocess
import pycs
import numpy as np
from .IO import *

__all__ = ['SLTimer', 'spl']

class SLTimer(object):
    '''
    Simple class for ingesting strong lens light curve data, and measuring the
    time .
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
        Download method to get the datafile from a url.
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
        Read in method to specify which type of datafile will be read in.
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
        return spl(self.lcs)

    #========================================================== Plotting light curves

    def display_light_curves(self,filename=None,jdrange=(None)):
        '''
        To display the lightcurves.
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
        self.lcs = whiten(self.lcs)
        return

    #===================================================== Microlensing

    def add_polynomial_microlensing(self):
        '''
        To add polynomial microlensing to each lightcurve.
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
        To add spline microlensing to each light curve.
        '''
        pycs.gen.splml.addtolc(self.lcs[0], knotstep=150)
        pycs.gen.splml.addtolc(self.lcs[1], knotstep=150)
        if self.Nim == 4:
            pycs.gen.splml.addtolc(self.lcs[2], knotstep=150)
            pycs.gen.splml.addtolc(self.lcs[3], knotstep=150)
        return

    #============================================= Primary workhorse method

    def estimate_time_delays(self,method='pycs',microlensing='spline',agn='spline',error=None):
        '''
        Provides both polynomial and spline time delays.
        '''
        if method == 'pycs':
            print "You are using the pycs method."
        else:
            print "The only available method is 'pycs' - exiting."
            return

        # Tell the lightcurves that their model is going to include microlensing:
        if microlensing == 'polynomial':
            self.add_polynomial_microlensing()
        elif microlensing == 'spline':
            self.add_spline_microlensing()
        else:
            pass
        # Keep a record:
        self.microlensing = microlensing

        # Optimize the model for both microlensing and intrinsic variability:
        if agn == 'spline':
            self.agn = self.optimize_spline_model()
        else:
            print "Error: only free-knot spline models are available for AGN variability at present."
            return

        # Print out time delays:
        self.report_time_delays()

        # Do error analysis, if required:
        if error == 'complete':
            self.estimate_uncertainties()
        elif error == 'intrinsic variance':
            self.find_intrinsic_variance()
        else:
            return
    
    #===================================================== Evaluate the fitting

    def computeLikelihood_MC(self,nsample=1000,nprocess=5,\
	rangeList=[[-100,100],[-100,100],[-100,100]]):
       	'''
        compute the likelihood by Montecarlo method
        '''
        import corner
        from multiprocessing import Pool
        from functools import partial
	import matplotlib
	# Force matplotlib to not use any Xwindows backend.
	matplotlib.use('Agg')
        ndim = 3
        dAB=np.random.uniform(rangeList[0][0],rangeList[0][1],nsample)
        dAC=np.random.uniform(rangeList[1][0],rangeList[1][1],nsample)
        dAD=np.random.uniform(rangeList[2][0],rangeList[2][1],nsample)
        sample=np.vstack((dAB,dAC,dAD)).T
        p = Pool(processes=nprocess)
        chisquare=np.array(p.map(partial(getChiSquare,self.lcs),sample))
        weight=np.exp((chisquare-np.max(chisquare)))
	weight/=np.sum(weight)
        print("weighted time delays (dAB,dAC,dAD)(days) :",weight.T.dot(sample))
        fig=corner.corner(sample,labels=[r'$\Delta t_{AB}(days)$',r'$\Delta t_{AC}(days)$',r'$\Delta t_{AD}(days)$'],
                        weights=weight,plot_contours=True,
                        plot_density=True,hist_kwargs={"log":True})
        fig.savefig("likelihood_%s_samples.png"%nsample)
	

    #===================================================== Resimulating the Data

    def delete_old_files(self):
        '''
        To delete the old files from prior runs through the data.
        '''
        subprocess.call('rm -rfv sims_copies sims_mocks', shell=True)
        subprocess.call('rm -rfv sims_copies_opt_spl sims_copies_opt_disp sims_copies_opt_regdiff', shell=True)
        subprocess.call('rm -rfv sims_mocks_opt_spl sims_mocks_opt_disp sims_mocks_opt_regdiff', shell=True)
        print "The old files have been deleted."
        return

    def make_plain_copies(self,n=None,npkl=None):
        '''
        To make copies of the data.
        '''
        Ncopies = n*npkl
        print "Making",Ncopies,"copies of the original dataset:"
        pycs.sim.draw.multidraw(self.lcs, onlycopy=True, n=n, npkl=npkl, simset="copies")
        return

    def make_mock_light_curves(self,n=None,npkl=None):
        # (modellcs, modelspline) = pycs.gen.util.readpickle("optspline.pkl")
        '''
        To make mock lightcurves to provide an basis for the estimate uncertainties.
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
        print "Making",Nmocks,"synthetic datasets, varying time delays by +/-",truetsr/2.0,"days"
        pycs.sim.draw.saveresiduals(modellcs, modelspline)
        pycs.sim.draw.multidraw(modellcs, modelspline, n=n, npkl=npkl, simset="mocks",
                truetsr=truetsr, tweakml=[Atweakml, Btweakml, Ctweakml, Dtweakml])
        return

    #===================================================== Making Multiple Model Fits

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

    #========================================================= Error Analysis

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

    def estimate_uncertainties(self,n=None,npkl=None):
        self.delete_old_files()
        self.make_plain_copies(n=n,npkl=npkl)
        self.make_mock_light_curves(n=n,npkl=npkl)
        # Add in an option to use regdiff and disp here
        self.make_spline_model_fits_of_plain_copies()
        self.make_spline_model_fits_of_mock_light_curves()
        self.plot_intrinsic_variance_histograms()
        self.error_summary()
        return

    def find_intrinsic_variance(self,n=None,npkl=None):
        self.make_plain_copies(n=n,npkl=npkl)
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

def splNoShift(lcs):
    spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=50,verbose=False,shifttime=False)
    spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=30,verbose=False,shifttime=False)
    spline = pycs.spl.topopt.opt_fine(lcs, nit=10, knotstep=20,verbose=False,shifttime=False)
    return spline
####To compute the chisquare
def getChiSquare(lcs_original,delay):
    import copy
    lcs=copy.deepcopy(lcs_original)
    for l in lcs:
        l.resetshifts()
        l.resetml()
    for index, l in enumerate(lcs):
        pycs.gen.splml.addtolc(l,knotstep=150)
        if index!=0:
           l.timeshift=delay[index-1]
    spline = splNoShift(lcs)
    return spline.lastr2nostab 


