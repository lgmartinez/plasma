import numba
from numba import jit
#
@jit(nopython=True)
def get_report_interval(report_interval1,report_interval2,itermax,step_it,step_it2):
    '''
    defines the reporting intervals for data output in the 
    solution module

    Parameters:
    -----------
    itermax: maximum number of iterations for which to report on
    step_it: reporting step, on-screen summary
    step_it2: reporting step, print graphs
    '''
    report_interval_size = int((itermax)/step_it)
    for j in range(report_interval_size):
        report_interval1[j] = j*step_it
        report_interval2[j] = j*step_it2
    #---
    return report_interval1,report_interval2
