#this script calculates compressibility of plasma
#Luis Martinez
#October 2, 2016
import numpy
def get_compressibility(pressure,reference_pressure,density,reference_density):
    ''' 
    Defines compressibility
    Beta = (1/rho)*drho/dP

    Parameters:
    -----------
    pressure
    reference pressure
    density
    reference density

    Returns:
    --------
    Beta: plasma compressibility
    '''
    Beta = numpy.zeros_like(density)
    px = pressure.copy()
    pxref = reference_pressure.copy()
    rhox = density.copy()
    rhoxref = reference_density.copy()
    #
    for i in range(Beta.size):
        if ( (px[i]-pxref[i]) < 1e-20 ):
            Beta[i] = 0.
        else:
            Beta[i] = (1/rhox[i])*(rhox[i]-rhoxref[i])/(px[i]-pxref[i])
    #
    return Beta
