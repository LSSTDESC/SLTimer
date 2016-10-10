import pycs
import os, sys
import numpy as np
from pycs.gen.lc import factory

__all__ = ['read_in_rdb_data', 'read_in_tdc2_data', 'tdc2import', 'factory', 'flexibleimport', 'flux2magnitude', 'whiten', 'SilentOperation']

# ======================================================================

def read_in_rdb_data(datafile):
    """
    Reads in the datafiles in the rdb format.
    This format is used in the original PyCS tutorial.
    """
    lcs =  [
            pycs.gen.lc.rdbimport(datafile, 'A', 'mag_A', 'magerr_A', "Trial"),
            pycs.gen.lc.rdbimport(datafile, 'B', 'mag_B', 'magerr_B', "Trial"),
            pycs.gen.lc.rdbimport(datafile, 'C', 'mag_C', 'magerr_C', "Trial"),
            pycs.gen.lc.rdbimport(datafile, 'D', 'mag_D', 'magerr_D', "Trial"),
            ]
    return lcs

# ======================================================================

def read_in_tdc2_data(datafile,whiten=False):
    """
    Reads in the datafiles that will be used for the Time Delay
    Challenge 2. Has an option to take the data from multiple filters
    and "whiten" it to look like it came from a single filter.
    """

    Nim = count_images(datafile)

    lcs = [
            tdc2import(datafile, 'A', 'flux_A', 'flux_A_err',
                       "Image", units='nmgy'),
            tdc2import(datafile, 'B', 'flux_B', 'flux_B_err',
                       "Image", units='nmgy')
          ]

    if Nim == 4:
        lcs.append( tdc2import(datafile, 'C', 'flux_C', 'flux_C_err',
                               "Image", units='nmgy') )
        lcs.append( tdc2import(datafile, 'D', 'flux_D', 'flux_D_err',
                               "Image", units='nmgy') )

    if whiten:
        whiten(lcs)

    return lcs

# ===================================================================

def flexibleimport(filepath, jdcol=1, magcol=2, errcol=3, startline=8, flagcol=None, propertycols=None, telescopename="Unknown", object="Unknown", plotcolour="red", verbose = True, units='ABmags'):
    """
    The general form of file reading. We read only one lightcurve object.
    To comment a line in the input file, start the line with # like in python.

    :param jdcol: The column number containing the MHJDs. First column is number 1, not 0 !
    :param magcol: The column number containing the magnitudes.
    :param errcol: The column number containing the magnitude errorbars.
    :param flagcol: (default : None) The colum containing a mask (if available). This one should contain False (= will be masked) or True (not masked).

    :param propertycols: (default : None) is a dict : ``propertycols = {"fwhm":7, "skylevel":8}`` means that col 7 should be read in as property "fwhm" and 8 as "skylevel".
    :type propertycols: dictionary
    """

    if verbose : print "Reading \"%s\"..." % (os.path.basename(filepath))
    datafile = open(filepath, "r")
    filelines = datafile.readlines()[startline-1:] # we directly "skip" the first lines of eventual headers ...
    datafile.close()
    jds = []
    mags = []
    magerrs = []
    flags = []
    properties = []
    # Check all the lines for headers, empty lines, and consistent numbers of columns:

    elsperline = None
    for i, line in enumerate(filelines) :

        # Check for header lines:
        if line[0] == "#":
            print "Skipping header line %i : %s" % (i+startline, repr(line))
            continue

        # Check for "empty" lines - this does not seem to work correctly...
        if len(line.split()) < 5:
            print "Skipping empty line %i : %s" % (i+startline, repr(line))
            continue

        # Check the consistency of the number of columns:
        elements = line.split() # line is a string, elements is a list of strings

        if elsperline is not None:
            if len(elements) is not elsperline:
                raise RuntimeError, "Parsing error in line %i, check columns : \n%s" % (i+startline, repr(line))

        elsperline = len(elements)

        # Now, extend data arrays using elements of this row:
        jds.append(float(elements[jdcol-1]))
        x,xerr = float(elements[magcol-1]), float(elements[errcol-1])
        if units=='ABmags':
            mags.append(x)
            magerrs.append(xerr)
        elif units=='nmgy':
            m,merr = flux2magnitude(x,xerr)
            mags.append(m)
            magerrs.append(merr)
        else:
            raise RuntimeError, "Unknown units '%s'\n" % (units)

        # Also append flag array, as required:
        if flagcol is not None :
            strflag = str(elements[flagcol-1])
            if strflag == "True":
                flags.append(True)
            elif strflag == "False":
                flags.append(False)
            else:
                print "Flag error in line %i : %s" % (i+startline, repr(line))
        else:
            flags.append(True)

        propdict = {} # an empty dict for now
        if propertycols is not None :
            for (propname, propcol) in propertycols.items():
                propdict[propname] = str(elements[propcol-1])
        properties.append(propdict)

    #Call the factory function from pycs
    pycs.gen.lc.factory(jds, mags)

	# Make a brand new lightcurve object:
    newlc = factory(np.array(jds), np.array(mags), magerrs=np.array(magerrs), telescopename=telescopename, object=object)
    newlc.properties = properties[:]
    newlc.mask = np.array(flags[:])
    newlc.plotcolour = plotcolour
    nbmask = np.sum(newlc.mask==False)
    commentcols = "(%i, %i, %i)" % (jdcol, magcol, errcol)
    newlc.commentlist.insert(0, "Imported from %s, columns %s" % (os.path.basename(filepath), commentcols))
    if verbose: print "%s with %i points imported (%i of them masked)." % (str(newlc), len(newlc.jds), nbmask)
    return newlc

# ======================================================================

def flux2magnitude(x,xerr):
    """
    Only used if the units are nmgy
    Arguments x and xerr are inputs from the original datafile
    Make sure to convert flux to magnitude before whitening to increase precision
    """
    m = 22.5 - 2.5*np.log10(x)
    x_lower, x_upper = x - xerr, x + xerr
    if x_lower < 0 or x < 0:
        m=99.0
        merr=99.0
    else:
        m_upper = 22.5 - 2.5*np.log10(x_lower)
        m_lower = 22.5 - 2.5*np.log10(x_upper)
        merr = (m_upper - m_lower)/2
    return m,merr

# ======================================================================

def count_images(filename):
    """
    Determine whether a TDC2 datafile contains 2 or 4 images.
    """
    tdc2file = open(filename, "r")
    lines = tdc2file.readlines()
    tdc2file.close()
    Ncols = len(lines[-1].split())
    Nim = (Ncols - 2) / 2
    if Nim == 2 or Nim == 4:
        pass
    else:
        raise ValueError("Unexpected number of images ",Nim)
    return Nim

# ======================================================================

def tdc2import(filepath, object="Unknown", magcolname="flux",
               magerrcolname="flux_err", telescopename="Unknown",
               plotcolour="red", mhjdcolname="MJD", flagcolname = None,
               propertycolnames = ["band"], verbose = True,
               units='ABmags'):
    """
    A way to import lightcurves using tdc2 data files. Column indices
    and delimeters do not matter, just provide the column headers that you want to read.

    Propertycolnames is a list of column names that can be added as
    properties.
    Possible settings :
    None : Will not import any properties
    ["property1", "property2"] : Provide a list of properties to import.

    The default column names are "MJD", "flux", "flux_err", "mask".

    Flexibleimport is used under the hood.
    """
    if verbose : print "Checking header of \"%s\"..." % (os.path.basename(filepath))
    tdc2file = open(filepath, "r")
    tdc2filelines = tdc2file.readlines()
    tdc2file.close()
    # Loop over lines until we reach the last one that starts with a '#'
    firstchar, k = '#', 0
    while firstchar == '#':
        firstchar = tdc2filelines[k][0]
        k += 1
    headers = tdc2filelines[k-2].split()
    # print "headers: ", headers
    # e.g.
    # headers:  ['#', 'MJD', 'band', 'flux_A', 'flux_A_err', 'flux_B', 'flux_B_err']
    # Note the hash mark: this makes the MJD column number 1 already

    # We check if the headers you want are available :
    checknames = [mhjdcolname, magcolname, magerrcolname]
    if flagcolname : checknames.append(flagcolname)
    if propertycolnames : checknames.extend(propertycolnames)

    for name in checknames:
        if name not in headers:
            raise RuntimeError, 'I cannot find a column named "%s" in your file !' % (name)

    jdcol = headers.index(mhjdcolname) # No need to +1
    magcol = headers.index(magcolname) # No need to +1
    errcol = headers.index(magerrcolname) # No need to +1

    if flagcolname is not None :
        flagcol = headers.index(flagcolname) + 1
    else:
        flagcol = None

    if propertycolnames is not None :
        propertycols = {}
        for propertycolname in propertycolnames:
            propertycols[propertycolname] = headers.index(propertycolname) + 1
    else:
        propertycols = None

    newlc = flexibleimport(filepath, jdcol=jdcol, magcol=magcol, errcol=errcol, startline=20, flagcol=flagcol, propertycols={"band":2}, telescopename=telescopename, object=object, verbose=verbose, units=units)
    newlc.plotcolour = plotcolour
    return newlc

# ======================================================================

def mean_and_scatter(lcs):
    """
    Simple function to measure mean and standard deviation of eachlight
    curve in the list.
    """
    Nim = len(lcs)
    names = ['A', 'B', 'C', 'D']
    mean, scatter = dict(), dict()
    for j in range(Nim):
        mean[names[j]] = np.mean(lcs[j].mags)
        scatter[names[j]] = np.std(lcs[j].mags)
    return mean,scatter

# ======================================================================

def whiten(lcs):
    """
    Function to whiten the lightcurve data in the case of multiple
    filters.
    """
    mu, sigma = mean_and_scatter(lcs)
    print "whiten: before whitening, means =", mu
    print "whiten: before whitening, scatters =", sigma

    # How many different images do we have?
    Nim = len(lcs)

    # Get list of observation filters, and the names of the bands:
    filters = np.array([])
    Nobs = len(lcs[0].properties)
    for k in range(Nobs):
        filters = np.append(filters,lcs[0].properties[k]['band'])
    bands = np.unique(filters)
    print "whiten: detected bands:",str(bands)

    # Which observations were taken in which band?
    index = dict()
    for band in bands:
        index[band] = np.where(filters == band)[0]

    # Loop over images:
    for j in range(Nim):

        # Find mean mag in each filter's light curve:
        indiv_mean = dict()
        for band in bands:
            indiv_mean[band] = np.mean(lcs[j].mags[index[band]])

        # Now find the overall mean, and compute each filter's offsets:
        overall_mean = np.mean(lcs[j].mags)

        offset = dict()
        for band in bands:
            offset[band] = indiv_mean[band] - overall_mean

        # Finally, apply the offsets:
        for band in bands:
            lcs[j].mags[index[band]] -= offset[band]

    mu, sigma = mean_and_scatter(lcs)
    print "whiten: after whitening, means =", mu
    print "whiten: after whitening, scatters =", sigma

    return lcs

# ======================================================================

class SilentOperation(object):
    """
    Redirects stdout and stderr to /dev/null, if required.

    Parameters
    ----------
    stdout: stream
         Destination of standard output
    stderr: stream
         Destination of standard errors

    Notes
    -----
    This code was cribbed from @tstone2077 at http://codereview.stackexchange.com/questions/25417/is-there-a-better-way-to-make-a-function-silent-on-need
    """
    def __init__(self, stdout=None, stderr=None):
        self.devnull = open(os.devnull,'w')
        self._stdout = stdout or self.devnull or sys.stdout
        self._stderr = stderr or self.devnull or sys.stderr

    def __enter__(self):
        self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
        self.old_stdout.flush(); self.old_stderr.flush()
        sys.stdout, sys.stderr = self._stdout, self._stderr

    def __exit__(self, exc_type, exc_value, traceback):
        self._stdout.flush(); self._stderr.flush()
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr
        self.devnull.close()
