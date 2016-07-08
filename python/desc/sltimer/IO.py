import pycs
import os
import numpy as np
from pycs.gen.lc import factory

__all__ = ['read_in_rdb_data', 'read_in_tdc2_data', 'tdc2import', 'factory', 'flexibleimport','flux2magnitude']

def read_in_rdb_data(datafile):
    lcs =  [
            pycs.gen.lc.rdbimport(datafile, 'A', 'mag_A', 'magerr_A', "Trial"),
            pycs.gen.lc.rdbimport(datafile, 'B', 'mag_B', 'magerr_B', "Trial"),
            pycs.gen.lc.rdbimport(datafile, 'C', 'mag_C', 'magerr_C', "Trial"),
            pycs.gen.lc.rdbimport(datafile, 'D', 'mag_D', 'magerr_D', "Trial"),
            ]
    return lcs

def read_in_tdc2_data(datafile):
    lcs = [
            tdc2import(datafile, 'A', 'flux_A', 'flux_A_err', "Image", units='nmgy'),
            tdc2import(datafile, 'B', 'flux_B', 'flux_B_err', "Image", units='nmgy'),
            tdc2import(datafile, 'C', 'flux_C', 'flux_C_err', "Image", units='nmgy'),
            tdc2import(datafile, 'D', 'flux_D', 'flux_D_err', "Image", units='nmgy'),
            ]
    whiten(lcs)
    return lcs

#===================================================================unfinished

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
            print("Skipping header line %i : %s" % (i+startline, repr(line)))
            continue

        # Check for "empty" lines - this does not seem to work correctly...
        if len(line.split()) < 5:
            print("Skipping empty line %i : %s" % (i+startline, repr(line)))
            continue

        # Check the consistency of the number of columns:
        elements = line.split() # line is a string, elements is a list of strings

        if elsperline != None:
            if len(elements) != elsperline:
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

def tdc2import(filepath, object="Unknown", magcolname="flux", magerrcolname="flux_err", telescopename="Unknown", plotcolour="red", mhjdcolname="MJD", flagcolname = None, propertycolnames = ["band"], verbose = True, units='ABmags'):
    """
    A way to import lightcurves using tdc2 data files. Column indices and delimeters do not matter,
    just provide the column headers that you want to read.

    Propertycolnames is a list of column names that can be added as properties.
    Possible settings :
    None : Will not import any properties
    ["property1", "property2"] : Provide a list of properties to import.

    The default column names are "mhjd", "mag", "magerr", "mask".

    Flexibleimport is used under the hood.
    """
    if verbose : print "Checking header of \"%s\"..." % (os.path.basename(filepath))
    rdbfile = open(filepath, "r")
    rdbfilelines = rdbfile.readlines()
    rdbfile.close()
    headers = rdbfilelines[18].split() #--> code has to pass the headers and extract the headers; read the values into a dictionary
    #underlines = rdbfilelines[1].split()
    #if map(len, headers) != map(len, underlines): --> find another test
    #raise RuntimeError, "Error in parsing headers"
    #headerindices = np.array(range(len(headers))) + 1 # +1 as we use human convention in flexibleimport

    # We build the default property names :
    #if propertycolnames == "lcmanip": # Then we put it the default set, but only if available.
        #propertycolnames = set(["telescope", "fwhm", "ellipticity", "airmass", "relskylevel", "normcoeff", "nbimg"]).intersection(set(headers))
    #if propertycolnames == ["band"]:

    # We check if the headers you want are available :
    checknames = [mhjdcolname, magcolname, magerrcolname]
    if flagcolname : checknames.append(flagcolname)
    if propertycolnames : checknames.extend(propertycolnames)

    for name in checknames:
        if name not in headers:
            raise RuntimeError, 'I cannot find a column named "%s" in your file !' % (name)

    jdcol = headers.index(mhjdcolname) + 1 # +1 as we use human convention in flexibleimport
    magcol = headers.index(magcolname) + 1
    errcol = headers.index(magerrcolname) + 1

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

"""
def whiten(lcs):
    #filters = np.array([])
    for item in lcs[0:4]:
        print ('item',woo)

        '''
        for k in range(100):
            filters = np.append(filters,lcs[].properties[k]['band'])
            print(filters)
    return
            index = np.where(filters == 'u')
            print(filters[index],lcs[0].mags[index])
            Image_A_u_band_mean = (np.mean(lcs[0].mags))
    #print(np.unique(filters))

    return lcs
"""
