# ======================================================================
# License info here?
# ======================================================================
from __future__ import absolute_import
import os
import urllib
import subprocess
import pycs
import numpy as np
from .reading import *
from matplotlib import pyplot as plt
import matplotlib
import scipy as sp
import sys
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

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
        self.ml_knotstep = 350
        self.knotstep = 20
        self.Hbar = 70.
        self.sigmaH = 7.
        self.phibar = None
        self.sigmaPhi = None
        self.Q = 0
        self.sigma_intrinsic = 0
        self.noise_rescaled=False
        return

    def rescale_noise(self):
        if self.noise_rescaled:
            print("you cannot rescale noise twice")
            return
        print("add additional noise {0}".format(self.sigma_intrinsic))
        for lc in self.lcs:
            lc.magerrs = np.sqrt(self.sigma_intrinsic**2 + lc.magerrs**2)
        self.noise_rescaled = True

    def reset_noise(self):
        if self.noise_rescaled is False:
            print("you cannot rest before rescale noise")
            return
        print("delete additional noise {0}".format(self.sigma_intrinsic))
        for lc in self.lcs:
            lc.magerrs = np.sqrt(-self.sigma_intrinsic**2 + lc.magerrs**2)
        self.noise_rescaled = False
            
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
            Q_FP_ERR = get_tdc2_header(self.datafile)
            self.Q = Q_FP_ERR['Q']
            self.phibar = Q_FP_ERR['FP']
            self.sigmaPhi = Q_FP_ERR['FPErr']
        else:
            raise ValueError('Unrecognized or null format '+str(format))

        self.Nim = len(self.lcs)

        return

    def prior(self, t, positive_H=False):
        t=-t ##Because the time convention is different in PyCS and TDC2
        Hbar = self.Hbar
        sigmaH = self.sigmaH
        phibar = self.phibar
        sigmaPhi = self.sigmaPhi
        Q = self.Q/(3.0*1E5)
   # print(Q*phibar/Hbar)
        f = 1./(2*sigmaH*sigmaPhi*np.pi*Q)
        s = -(Hbar)**2/(sigmaH**2)+(-phibar**2)/(sigmaPhi**2)
        t = ((Hbar/(sigmaH**2)+(phibar*t)/(Q*sigmaPhi**2))**2)/(1./(sigmaH**2)+(t**2)/((sigmaPhi**2)*(Q**2)))
        normalize = np.max(t)+s
        m = np.exp(s+t-normalize)
        b = (Hbar/sigmaH**2+(phibar*t)/(Q*(sigmaPhi**2)))/(1./sigmaH**2+t**2/((sigmaPhi**2)*(Q**2)))
        a = 1./sigmaH**2+t**2/((sigmaPhi**2)*(Q**2))
        if positive_H:
            ft = (np.exp(-a*b**2)+np.sqrt(np.pi*a)*b*(sp.special.erf(np.sqrt(a)*b)+1))/(2*a)
        else:
            ft = b*np.sqrt(np.pi/a)
        return f*m*ft
    def chisquare_to_loglikelihood(self, chisquare):
        number_of_data = 0
        lognoise_sum = 0
        for lc in self.lcs:
            number_of_data += len(lc)
            lognoise_sum += np.sum(np.log(lc.magerrs))
        print(lognoise_sum)
        return -1./2.*chisquare-number_of_data/2.*np.log(2*np.pi)-lognoise_sum

    def add_prior_to_sample(self, result):
        prior = self.prior(result['dt_AB'], positive_H=True)
        original = np.zeros((result['dt_AB'].shape[0], 2))
        log_prior = np.zeros((result['dt_AB'].shape[0], 2))
        combined = np.zeros((result['dt_AB'].shape[0], 2))

        original[:, 0] = result['dt_AB']
        log_prior[:, 0] = result['dt_AB']
        combined[:, 0] = result['dt_AB']
        if 'log_likelihood' not in result.keys():
            original[:, 1] = self.chisquare_to_loglikelihood(result['chisquare'])
        else:
            original[:, 1] = result['log_likelihood']
        log_prior[:, 1] = np.log(prior)
        combined[:, 1] = original[:, 1]+log_prior[:, 1]
        return [original, log_prior, combined]

    def optimize_spline_model(self):
        '''
        Optimizes a spline model for the intrinsic variability.
        '''
        return spl(self.lcs, knotstep=self.knotstep)

    #========================================================== Plotting light curves
    def display_light_curves(self, filename=None, jdrange=(None), title=None,
                             given_curve=None):
        '''
        Displays the lightcurves in a single panel plot.
        '''
        if given_curve is not None:
            if len(given_curve) == 2:
                lcs, agn = given_curve
            else:
                lcs = given_curve
                agn = None
        else:
            lcs = self.lcs
            agn = None
        pycs.gen.mrg.colourise(lcs)
        # Replace the following with an optional input list of shifts
        #lcs[1].shifttime(-5.0)
        #lcs[2].shifttime(-20.0)
        #lcs[3].shifttime(-70.0)
        pycs.gen.lc.display(lcs, [agn], figsize=(20, 7),
                            jdrange=jdrange, title=title, nicefont=True)
        # lcs = pycs.gen.util
        # for l in lcs:
        #     l.resetshifts()
        if filename is not None:
            pycs.gen.lc.display(lcs, [agn], figsize=(20, 7),
                                jdrange=jdrange, title=title, nicefont=True,
                                filename=filename)
        return

    def select_bands(self, bands):
        '''
        select bands you want to keep

        Notes:
        ------

        .. warning:: this function will change the light curve in SLTimer
        '''
        self.lcs = select_bands(self.lcs, bands)

    def reset_lc(self):
        for l in self.lcs:
            l.resetshifts()
            l.resetml()
        return

    def whiten(self, seasonal=False):
        '''
        Whitens a set of multi-filter light curves to a single fictitious band.
        '''
        if seasonal:
            self.lcs = whiten_season(self.lcs)
        else:
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
        spline_microlensing(self.lcs, self.ml_knotstep)
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

    #===================================================== Evaluate the fitting
    def compute_chisq(self, delay, batch=False, getlcs=False):
        """
        return chisquare of spline fitting given time delay

        Parameters
        ----------
        delay : 1D array
            array contains time delays for each light curve. The convention is
            [dt_AB, dt_AC, dt_AD]

        batch : bool
            if batch==True, then delay can be a two dimensional array with each
                row contains a set of time delay sample.

        """
        if batch:
            chisquare = []
            for item in delay:
                chisquare.append(get_chi_squared(
                                 lcs_original=self.lcs,
                                 ml_knotstep=self.ml_knotstep, delay=item,
                                 getlcs=False, knotstep=self.knotstep
                                 ))
            return chisquare
        return get_chi_squared(lcs_original=self.lcs,
                               ml_knotstep=self.ml_knotstep,
                               getlcs=getlcs,
                               delay=delay, knotstep=self.knotstep)

    def generate_random_sample(self, rangeList, nsample):
        ndim = len(self.lcs)
        #Generate samples
        if rangeList is None:
            rangeList = [[-100, 100]]*(ndim-1)
        d = []
        for item in xrange(ndim-1):
            d.append(np.random.uniform(rangeList[item][0], rangeList[item][1],
                     nsample))
        sample = np.array(d).T
        return sample

    def write_out_to(self, result, outName):
        file_name = "{0}_delay_chi2_{1}_samples.txt".format(outName,
                                                            result.shape[0])
        print result.shape[1]
        if result.shape[1]==4:
            names = ["AB"]
        else:
            names = ["AB", "AC", "AD"]
        header = "Smaples time delay for simple montecarlo and their corresponding \
        chisquare. \n"
        for item in names:
            header += "   dt_"+item
        header += "   chisquare"
        header += "   log_likelihood"
        header += "   sigma_intrinsic"
        np.savetxt(file_name, result, header=header, comments="# ")
        return

    def plot_likelihood_from_file(self, file_name,
                                  chisquare=False, likelihood=False, bins=20,
                                  outName="from_file_", corner_plot=True,
                                  add_prior=False, batch_sigma=False,
                                  method="plot_log_file"):
        file = open(file_name)
        lines = file.readlines()
        file.close()
        keys = " ".join(lines[1].split()).split()[1:]
        result = {}
        for key in keys:
            result[key] = np.array([])
        for l in lines[2:]:
            array = l.split()
            for index, key in enumerate(keys):
                result[key] = np.append(result[key], eval(array[index]))
        if add_prior:
            result_new = self.add_prior_to_sample(result)
            if method == "plot exp in same graph":
                self.plot_explikelihood_same_file(result_new,
                                                  outName+file_name[-10:],
                                                  bins=bins)

            else:
                if batch_sigma:
                    result_new.append(result['sigma_intrinsic'])
                self.plot_log_likelihood_with_prior(result_new,
                                                    outName+file_name[-10:],
                                                    bins=bins,
                                                    batch_sigma=batch_sigma)
        else:
            result_new = []
            for keys in result.keys():
                if 'dt' in keys:
                    result_new.append(result[keys])
            if likelihood:
                result_new.append(result['log_likelihood'])
            else:
                result_new.append(result['chisquare'])
            result_new = np.array(result_new)
            result_new = result_new.T
            self.plot_likelihood(result_new, outName+file_name[-10:],
                                 chisquare=chisquare, bins=bins,
                                 corner_plot=corner_plot, likelihood=likelihood)
        return

    def plot_explikelihood_same_file(self, result, outName,
                                     bins=20):
        fig = plt.figure(figsize=(5, 10))
        ax = fig.add_subplot(111)
        for item in result:
            exp_prop = np.exp(item[:, -1]-np.max(item[:, -1]))
            item[:, -1] = exp_prop/(np.sum(exp_prop)*(item[1, 0]-item[0, 0]))
            print np.sum(exp_prop)
        original = result[0]
        prior = result[1]
        combined = result[2]
        labelcolor = ["likelihood", "b"]
        self.internal_plot(result=original,
                           bins=bins, corner_plot=False,
                           ax=ax, chisquare=True, labelcolor=labelcolor)
        labelcolor = ["prior", "k"]
        self.internal_plot(result=prior,
                           bins=bins, corner_plot=False,
                           ax=ax, chisquare=True, labelcolor=labelcolor)

        labelcolor = ["posterior", "r"]
        self.internal_plot(result=combined,
                           bins=bins, corner_plot=False,
                           ax=ax, chisquare=True, labelcolor=labelcolor)
        ax.set_ylabel(r"probability")
        plt.legend(loc=(1.05, 0.9))
        fig.suptitle("prior_likelihood_posterior")
        fig.savefig("{0}_prior_likelihood_posterior_{1}_samples.png".format(
                    outName, result[0].shape[0]))
        return



    def plot_log_likelihood_with_prior(self, result, outName,
                                       bins=20, batch_sigma=False):
        import matplotlib.gridspec as gridspec
        if batch_sigma:
            sigmaArr = np.unique(result[3])
            sigmas = result[3]
            originals = []
            priors = []
            combineds = []
            for sigma in sigmaArr:
                index = np.where(sigmas == sigma)
                originals.append(result[0][index])
                priors.append(result[1][index])
                combineds.append(result[2][index])
        fig = plt.figure(figsize=(5, 10))
        gs = gridspec.GridSpec(3, 1, height_ratios=[4, 4, 4])
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[2, 0])
        gs.update(left=0.01, right=0.99, hspace=0.3)
        def plotBatch(axes, result_inside, labelcolor=None):
            ax1 = axes[0]
            ax2 = axes[1]
            ax3 = axes[2]
            original = result_inside[0]
            prior = result_inside[1]
            combined = result_inside[2]
            self.internal_plot(result=original,
                               bins=bins, corner_plot=False,
                               ax=ax1, chisquare=True, labelcolor=labelcolor)
            ax1.set_ylabel(r'$log(L)$')
            ax1.set_xlabel('')
            self.internal_plot(result=prior,
                               bins=bins, corner_plot=False,
                               ax=ax2, chisquare=True, labelcolor=labelcolor)
            ax2.set_ylabel(r'$log(prior)$')
            ax2.set_xlabel('')
            self.internal_plot(result=combined,
                               bins=bins, corner_plot=False,
                               ax=ax3, chisquare=True, labelcolor=labelcolor)
            ax3.set_ylabel(r'$log(posterior)$')
        axes = [ax1, ax2, ax3]
        if batch_sigma:
            import matplotlib.cm as cm
            colors = iter(cm.rainbow(np.linspace(0, 1, len(sigmaArr))))
            for index, sigma in enumerate(sigmaArr):
                result_new = [originals[index], priors[index], combineds[index]]
                labelcolor = ["$\sigma_{int}=$"+str(sigma), next(colors)]
                plotBatch(axes, result_new, labelcolor=labelcolor)
        else:
            plotBatch(axes, result)
        plt.legend(loc=(1.05, 2.9))
        fig.suptitle("log likelihood")
        fig.savefig("{0}_likelihood_{1}_samples.png".format(outName,
                                                            result[0].shape[0]))
        return

    def plot_likelihood(self, result, outName, plot_contours=True,
                        plot_density=True, chisquare=False, bins=20,
                        corner_plot=True, ax=None, likelihood=False):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        newFig = self.internal_plot(result=result, plot_contours=plot_contours,
                                    plot_density=plot_density,
                                    chisquare=chisquare,
                                    bins=bins, corner_plot=corner_plot, ax=ax,
                                    likelihood=likelihood)
        if newFig is not None:
            fig = newFig
        if likelihood:
            title = r"$log(likelihood) plot$"
        else:
            title = r"$\chi^2 plot$"
        fig.suptitle(title)
        fig.savefig("{0}_likelihood_{1}_samples.png".format(outName,
                    result.shape[0]))

    def internal_plot(self, result, plot_contours=True,
                      plot_density=True, chisquare=False, likelihood=False, bins=20,
                      corner_plot=True, ax=None, labelcolor=None):
        import corner
        sample = result[:, :-1]
        if not chisquare:
            weight = chi2_to_weight(result[:, -1])
        else:
            weight = result[:, -1]
        if corner_plot:
            fig = corner.corner(sample, bins=bins,
                                labels=[r'$\Delta t_{AB}(days)$',
                                        r'$\Delta t_{AC}(days)$',
                                        r'$\Delta t_{AD}(days)$'],
                                weights=weight, plot_contours=plot_contours,
                                plot_density=plot_density,
                                max_n_ticks=10,
                                use_math_text=True
                                )
        else:
            if sample.shape[1] != 1:
                print("corner=False can only be true when there is only 1D sample")
            sample = sample.ravel()
            bins = np.linspace(sample.min(), sample.max(), bins)
            mask = np.where(weight != -np.inf)
            wd, b = np.histogram(sample[mask], bins=bins, weights=weight[mask])
            counts, b = np.histogram(sample[mask], bins=bins)
            bincentres = [(b[i]+b[i+1])/2. for i in range(len(b)-1)]
            ax_min = max(sample.min(), -100)
            ax_max = min(sample.max(), 100)
            ax.set_xticks(np.linspace(ax_min, ax_max, 21), 5)
            ax.set_xlim(sample.min(), sample.max())
            ax.set_xlabel(r'$\Delta t_{AB}(days)$')
            if likelihood:
                ax.set_ylabel(r'$log(L)$')
            else:
                ax.set_ylabel(r'$\chi^2$')
            if labelcolor is not None:
                color = labelcolor[1]
                label = labelcolor[0]
            else:
                color = 'k'
                label = None
            ax.step(bincentres, wd/counts, where='mid', color=color,
                    linestyle="-", label=label)
            import matplotlib.ticker as mtick
            ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
            fig = None
        return fig

    def compute_likelihood_simpleMC(self, nsample=1000, nprocess=5,
                                    rangeList=None, outName="",
                                    save_file=True, samples=None):
        '''
        Compute the likelihood by Monte Carlo method
        '''
        from multiprocessing import Pool
        from functools import partial
        import time
        if samples is not None:
            sample = samples
            nsample = len(sample)
        else:
            sample = self.generate_random_sample(rangeList=rangeList,
                                                 nsample=nsample)
        #calculate the chisquare
        start = time.time()
        p = Pool(processes=nprocess)
        chisquare = np.array(p.map(partial(
                                    get_chi_squared,
                                    lcs_original=self.lcs,
                                    ml_knotstep=self.ml_knotstep,
                                    getlcs=False,
                                    knotstep=self.knotstep),
                                   sample))
        end = time.time()
        print("Multiprocessing used {0} seconds.".format(end-start))
        weight = chi2_to_weight(chisquare)
        print("min chisquare,", np.min(chisquare))
        print("#"*20)
        print("weighted time delays (dAB,dAC,dAD)(days) :",
              weight.T.dot(sample))
        sigma_intric = [self.sigma_intrinsic]*len(sample)
        results = np.column_stack((sample,
                                   chisquare,
                                   self.chisquare_to_loglikelihood(chisquare),
                                   sigma_intric))

        if save_file:
            self.write_out_to(results, outName)
            self.plot_likelihood(results, outName)
        return sample[np.argmin(chisquare)]

    def degree_of_freedom(self):
        spline = pycs.spl.topopt.opt_rough(self.lcs, nit=1,
                                           knotstep=self.knotstep,
                                           verbose=False)
        num = len(spline.t)
        spline = pycs.spl.topopt.opt_rough(self.lcs, nit=1,
                                           knotstep=self.ml_knotstep,
                                           verbose=False)
        num_ml = len(spline.t)
        free_param = num*2+4+len(self.lcs)*(num_ml*2+4)+4
        nm_constraint = 0
        for l in self.lcs:
            nm_constraint += len(l)
        print("knotstep for intrinsic fluctuation is: {0}".format(self.knotstep))
        print("knotstep for micro lensing is: {0}".format(self.ml_knotstep))
        print("number of data points is: {0}".format(nm_constraint))
        dof = nm_constraint-free_param
        return {"dof" : dof, "# data" : nm_constraint}

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

        elif method == 'simpleMC':
            bestGuess = self.compute_likelihood_simpleMC(nsample=10,
                                                         nprocess=4,
                                                         save_file=False)
            dt = {'AB': bestGuess[0]}
            if self.Nim == 4:
                dt = {'AC': bestGuess[1]}
                dt = {'AD': bestGuess[2]}
        elif method == 'prior':
            dt = {'AB': -1.*self.Q*self.phibar/(3.0*1E5*self.Hbar)}
            print(dt)
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
        return self.time_delays

# ======================================================================
# End of the SLTimer class.
# ======================================================================


# Optimizer functions (could go in "optimize.py" instead?)
def spl(lcs, shifttime=True, verbose=True, knotstep=20):
    spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=5/2.*knotstep,
                                       shifttime=shifttime, verbose=verbose)

    spline = pycs.spl.topopt.opt_rough(lcs, nit=5, knotstep=3/2.*knotstep,
                                       shifttime=shifttime, verbose=verbose)

    spline = pycs.spl.topopt.opt_fine(lcs, nit=10, knotstep=knotstep,
                                      shifttime=shifttime, verbose=verbose)
    return spline


def spline_microlensing(lcs, ml_knotstep):
    if ml_knotstep is None:
        print("you didn't add any microlensing")
    else:
        for l in lcs:
            pycs.gen.splml.addtolc(l, knotstep=ml_knotstep)
    return
# To compute the chisquare


def get_chi_squared(delay, lcs_original, ml_knotstep, getlcs, knotstep=20):
    import copy
    lcs = copy.deepcopy(lcs_original)
    for l in lcs:
        l.resetshifts()
        l.resetml()
    spline_microlensing(lcs, ml_knotstep)
    for index, l in enumerate(lcs):
        if index != 0:
            l.timeshift = delay[index-1]
    spline = spl(lcs, verbose=False, shifttime=False, knotstep=knotstep)
    if getlcs:
        return [lcs, spline]
    else:
        return spline.lastr2nostab


def chi2_to_weight(chisquare):
    weight = np.exp(-0.5*(chisquare-np.min(chisquare)))
    weight /= np.sum(weight)
    return weight
