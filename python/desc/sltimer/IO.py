import pycs
import os
import numpy as np

__all__ = ['read_in_rdb_data', 'read_in_tdc2_data', 'tdc2import', 'factory', 'flexibleimport']

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
    return lcs
    
#===================================================================unfinished

def factory(jds, mags, magerrs=None, telescopename="Unknown", object="Unknown", verbose=False):
    """
    Returns a valid lightcurve object from the provided arrays.
    The numpy arrays jds and mags are mandatory. If you do not specify a third array containing the magerrs,
    we will calculate them "automatically" (all the same value), to avoid having 0.0 errors.
        
    @type	jds: 1D numpy array
    @param	jds: julian dates
    @type	mags: 1D numpy array
    @param	mags: magnitudes
    @type	magerrs: 1D numpy array
    @param	magerrs: optional magnitude errors
    
    @todo: improve it and use this in file importing functions
    
    """
    # Make a brand new lightcurve object :
    newlc = pycs.gen.lc.lightcurve()
            
    # Of couse we can/should check a lot of things, but let's be naive :
            
    newlc.jds = np.asarray(jds)
    newlc.mags = np.asarray(mags)
            
    if magerrs is None:
        newlc.magerrs = np.zeros(len(newlc.jds)) + 0.1
    else:
        newlc.magerrs = np.asarray(magerrs)
            
    if len(newlc.jds) != len(newlc.mags) or len(newlc.jds) != len(newlc.magerrs):
        raise RuntimeError, "lightcurve factory called with arrays of incoherent lengths"
            
    newlc.mask = newlc.magerrs >= 0.0	# This should be true for all !
        
    newlc.properties = [{}] * len(newlc.jds)

    newlc.telescopename = telescopename
    newlc.object = object
        
    newlc.setindexlabels()
    newlc.commentlist = []
        
    newlc.sort() # not sure if this is needed / should be there
        
    #newlc.validate()
        
    if verbose: print "New lightcurve %s with %i points" % (str(newlc), len(newlc.jds))
    return newlc


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
        if flagcol != None :
            strflag = str(elements[flagcol-1])
            if strflag == "True":
                flags.append(True)
            elif strflag == "False":
                flags.append(False)
            else:
                print "Flag error in line %i : %s" % (i+startline, repr(line))
        else:
            flags.append(True)

        # Finally, extend the properties array:
        propdict = {}
        if propertycols is not None :
            for (propname, propcol) in propertycols.items():
                propdict[propname] = str(elements[propcol-1])
        properties.append(propdict)

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
    m = 22.5 - 2.5*np.log10(x)
    x_lower, x_upper = x - xerr, x + xerr
    if x_lower < 0 or x < 0:
        m=99.0
        merr=99.0
    else:
        merr = 22.5 - 2.5*np.log10(xerr)
    return m,merr


def tdc2import(filepath, object="Unknown", magcolname="flux", magerrcolname="flux_err", telescopename="Unknown", plotcolour="red", mhjdcolname="MJD", flagcolname = None, propertycolnames = ["band"], verbose = True, units='ABmags'):
    """
    The relaxed way to import lightcurves, especially those from cosmouline, provided they come as rdb files.
    Don't care about column indices, just give the column headers that you want to read.
    
    Propertycolnames is a list of column names that I will add as properties.
    Possible settings :
    "lcmanip" : (default) I will import the standard cols from lcmanip / cosmouline, if those are available.
    None : I will not import any properties
    ["property1", "property2"] : just give a list of properties to import.
    
    The default column names are those used by the method :py:meth:`pycs.gen.lc.lightcurve.rdbexport` : "mhjd", "mag", "magerr", "mask".
    
    We use flexibleimport under the hood.
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
    if propertycolnames == "lcmanip": # Then we put it the default set, but only if available.
        propertycolnames = set(["telescope", "fwhm", "ellipticity", "airmass", "relskylevel", "normcoeff", "nbimg"]).intersection(set(headers))

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
    
    if flagcolname != None :
        flagcol = headers.index(flagcolname) + 1
    else:
        flagcol = None

    if propertycolnames != None :
        propertycols = {}
        for propertycolname in propertycolnames:
            propertycols[propertycolname] = headers.index(propertycolname) + 1
    else:
        propertycols = None

    newlc = flexibleimport(filepath, jdcol=jdcol, magcol=magcol, errcol=errcol, startline=20, flagcol=flagcol, propertycols=propertycols, telescopename=telescopename, object=object, verbose=verbose, units=units)
    newlc.plotcolour = plotcolour
    return newlc


"""
#convert flux into magnitude; before we return
#M = -2.5log10f (ab maggies)
#Mk = 99.0; Magerror = 99.0
#if f < 0
#read the fluxes; lightcurves.mags = -2.5*np.log10(lc.mags)
#need an algorithm to convert flux error to mag error; suggest an algorithm
#F+sigma(equation); = m+sigma; F-sigma = m-sigma; then plug that into the equation
#Finite differencing;
#-2sigmam = -2.5log(f-sigma)+2.5log(f+sigma); divide both sides by two; if we still get a negative number; just take the positive part
#flux_to_magnitude --> flux to magnitude; flux errors to mag errors
#As input get an array  flux inputs and flux errors; return as an array the mag/magerror; code: mag, mag error = mag, mag error lc.mag
#fluxes will be stored in an array; lc.mag, lc.magerror --> could act on the object; lightcurve = convert flux to mag
#--readin needs to be able to do this when the; readintdc2data will call fleximport; have it call tdc2import
Whitening:
Offset in magnitude; factor of flux; comes from a number of place: magnitudes measure faintness, expect g band pts to be below r band pts; there is a constant factor
Turn the multi-filter light curve into effectively a single light curve; w band (for whitened); get rid of the offset, Mean of the r-band; then mean of the g-band; that distance is the difference between g and r; divide by two; subtract the higher pts by that; then bring the bottom by that; this will be the w-band
To do this: function called whiten; take in a list of lightcurve objects and return a lst of whitened lightcurve objects; overwrite lcs objects to whitened versions; set a flag in the timer, called self.whitened=True (as a reminder that you have whitened) be able to whitened individually of each

note:first convert to magnitudes then whiten
If you make the whitening funciton a null opt --> should get back the synthetic w band magnitudes; should see a bunch of offsets; light curves will be defined
lcs object will need a column for filter
pycs.lc.lightcurve(to import lightcurve)
---
preserve the pycs structure; my own version of factory and flexibleimport; flexibleimport --> multifilter=True; units=flux(instead of mags)
"""