# ==============================================================================
# License info here?
# ==============================================================================

import sltimer

# ==============================================================================

class SLTimer(object):
    '''
    Simple class for ingesting strong lens light curve data, and measuring the
    time delays.
    '''
    def __init__(self):
        return

    def read(self,lcfile):
        print("No IO enabled yet.")
        return

    def run(self,algorithm=None):
        if algorithm is None:
            print("No algorithms coded yet.")
        return

# ==============================================================================
# Need better tests...

if __name__ == '__main__':

    timer = sltimer.SLTimer()

    timer.read("lightcurve.txt")

    timer.run()

# ==============================================================================
