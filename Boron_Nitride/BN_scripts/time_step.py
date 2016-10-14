#this script provides the algorithm for the adaptive time step scheme
#Luis Martinez
#October 3, 2016
import numpy
def get_adaptive_time_step(cflenergy,cflmomentum,energyflux,momentumflux,\
                           delx,delt_ref,maximum_velocity):
    ''' Adaptive time step
    '''
    dtlim = 4e-9
    cfle_max = 125.
    cflm_max = 0.75
    #
    if (cflenergy>cfle_max) or (cflmomentum>cflm_max): 
        dte = (delx*cfle_max)/energyflux
        dtm = (delx*cflm_max)/momentumflux
        dtlim = numpy.minimum(dte,dtm)
    #
    if (dtlim>delt_ref):
        dte = (delx*cfle_max)/energyflux
        dtm = (delx*cflm_max)/momentumflux
        dtlim = numpy.minimum(dte,dtm)
    #
    return dtlim
