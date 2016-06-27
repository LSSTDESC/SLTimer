
def read_in_rdb_data(datafile):
    import pycs
    lcs =  [
            pycs.gen.lc.rdbimport(datafile, 'A', 'mag_A', 'magerr_A', "Trial"),
            pycs.gen.lc.rdbimport(datafile, 'B', 'mag_B', 'magerr_B', "Trial"),
            pycs.gen.lc.rdbimport(datafile, 'C', 'mag_C', 'magerr_C', "Trial"),
            pycs.gen.lc.rdbimport(datafile, 'D', 'mag_D', 'magerr_D', "Trial"),
           ]
    return lcs

def read_in_tdc2_data(datafile):
    import pycs
    lcs = [
            pycs.gen.lc.rdbimport(datafile, 'A', 'flux_A', 'flux_A_err', "Trial"),
            pycs.gen.lc.rdbimport(datafile, 'B', 'flux_B', 'flux_B_err', "Trial"),
            pycs.gen.lc.rdbimport(datafile, 'C', 'flux_C', 'flux_C_err', "Trial"),
            pycs.gen.lc.rdbimport(datafile, 'D', 'flux_D', 'flux_D_err', "Trial"),
            ]
    return lcs
