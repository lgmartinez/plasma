#Luis Martinez
#September 29, 2016
import numpy
from parameters import  k_B,k_N,k_Ni,k_Co,mB,mNi,mCo,mN
'''
This script contains the following functions:
thermophysical_properties.get_thermophysical_properties
 
This function obtains:
1. thermal conductivity
2. shear viscosity
3. specific heat
4. specific gas constant
'''
#
def get_thermophysical_properties(temperature,concentration1,concentration2,\
                                 concentration3,concentration4):
    '''
    finds thermal conductivity,viscosity,specific heat,and
    specific gas constant
    
    Parameters:
    -----------
    temperature: array of domain temperature
    concentration1: Nitrogen
    concentration2: Boron
    concentration3: Nickel
    concentration4: Cobalt

    Returns:
    --------
    mux: viscosity [Pa-sec]
    kx: thermal conductivity [W/m-K]
    cpx: speficic heat [J/Kg-K]
    Rsx: specific gas constant in domain [J/kg-K]
    '''
    #concentrations
    cs1 = concentration1.copy() #nitrogen
    cs2 = concentration2.copy() #boron
    cs3 = concentration3.copy() #nickel
    cs4 = concentration4.copy() #cobalt
    #viscosity
    mux = numpy.zeros_like(temperature) 
    mus1 = numpy.zeros_like(temperature)
    mus2 = numpy.zeros_like(temperature)
    mus3 = numpy.zeros_like(temperature)
    mus4 = numpy.zeros_like(temperature) 
    #thermal conductivity
    kx = numpy.zeros_like(temperature) 
    ks1 = numpy.zeros_like(temperature) 
    ks2 = numpy.zeros_like(temperature)
    ks3 = numpy.zeros_like(temperature) 
    ks4 = numpy.zeros_like(temperature) 
    #specific heat
    cpx = numpy.zeros_like(temperature)
    cps1 = numpy.zeros_like(temperature)
    cps2 = numpy.zeros_like(temperature) 
    cps3 = numpy.zeros_like(temperature) 
    cps4 = numpy.zeros_like(temperature) 
    #
    Tx = temperature.copy() # temperature
    array_size = numpy.size(temperature)
    #
    #--------thermal conductivity---------
    #conductivity of Nitrogen
    aN4,bN4 = 1.1973364764750705e-11,-1.0022238377799794e-07 
    cN4,dN4,eN4 = 0.00025655111953739039,-0.1423732005907308,33.464978104210182
    ks1[:] = 1e-3*(aN4*(Tx[:]**4) + bN4*(Tx[:]**3) +\
                   cN4*(Tx[:]**2) + dN4*Tx[:] + eN4)
    #conductivity of Boron
    ks2[:] = k_B
    #conductivity of Nickel
    ks3[:] = k_Ni
    #conductivity of Cobalt
    ks4[:] = k_Co
    # total viscosity
    kx[:] = cs1[:]*ks1[:] + cs2[:]*ks2[:] + cs3[:]*ks3[:] + cs4[:]*ks4[:] 
    # kx[-1] = kth[-2]
    #--------dynamic viscosity------------------
    #viscosity of nitrogen
    aNmu,bNmu = 4.2177299011972125e-09, -4.8160156609647342e-05
    cNmu,dNmu = 0.40073205811284623, 44.208627576799387
    mus1[:] = 1e-7*(aNmu*(Tx[:]**3) + bNmu*(Tx[:]**2) + cNmu*(Tx[:]) + dNmu)
    #viscosity of boron (gas)
    mus2[:] = 5.4e-5 #viscosity used in carbon-helium simulation #1e-5 #2e-3
    #viscosity of nickel (gas)
    mus3[:] = 0.
    #viscosity of cobalt (gas)
    mus4[:] = 0.
    #total viscosity
    mux[:] = cs1[:]*mus1[:] + cs2[:]*mus2[:] + cs3[:]*mus3[:] + cs4[:]*mus4[:]
    #mu[-1] = mu[-2]#2*mu[-2] - mu[-3]
    #---------specific heat ---------------
    #--- nitrogen
    for i in range(array_size):
        cps1[i] = 6.66 + 1.02*(1e-3)*Tx[i] # {cal}/{gram-mole K}
    #conversion:
    mw_N2 = 2*14.0067*(1e-3) # molar mass in kg
    cps1[:] = cps1[:] * (4.184)/mw_N2
    #--- Boron
    for i in range(array_size):
        if (Tx[i] < 1200.):
            cps2[i] = 1.54 + 4.49*(1e-3)*Tx[i]
        else:
            cps2[i] = 1.54 + 4.49*(1e-3)*1200.
    #conversion
    mw_B = 10.811*1e-3 #kg
    cps2[:] = cps2[:] * 4.814/mw_B
    #--- nickel
    for i in range(array_size):
        if ( Tx[i] < 633.0):
            cps3[i] = 4.06 + 7.04*(1e-3)*Tx[i]
        elif (Tx[i] >= 633. and Tx[i] < 1725.):
            cps3[i] = 6.0 + 1.8*(1e-3)*Tx[i]
        elif (Tx[i] >=1725):
            cps3[i] = 9.20
    #conversion
    mw_Ni = 58.6934*1e-3
    cps3[:] = cps3[:] * (4.184)/mw_Ni
    #--- cobalt
    for i in range (array_size):
        if (Tx[i] < 718.):
            cps4[i] = 4.72 + 4.3*(1e-3)*Tx[i]
        elif ( Tx[i] >= 718. and  Tx[i] < 1400.):
            cps4[i] = 3.30 + 5.86*(1e-3)*Tx[i]
        elif (Tx[i] >= 1400.):
            cps4[i] = 9.60
    #conversion:
    mw_Co = 58.933195*1e-3
    cps4[:] = cps4[:] * 4.184/mw_Co
    # total specific heat
    cpx[:] = cs1[:]*cps1[:] + cs2[:]*cps2[:] +\
             cs3[:]*cps3[:] + cs4[:]*cps4[:]
    #
    #----------specific gas constant--------------
    Ru = 8134.
    Rsx = numpy.zeros_like(concentration1)
    molarmass = numpy.zeros_like(concentration1)
    #molar mass
    ms1 = mN
    ms2 = mB
    ms3 = mNi
    ms4 = mCo
    #total molar mass
    molarmass[:] = 1./((cs1[:]/ms1) + (cs2[:]/ms2) +\
                       (cs3[:]/ms3) + (cs4[:]/ms4))
    #gas constant
    Rsx[:] = Ru/molarmass[:]
    #----------------------------------------------
    #
    return kx,mux,cpx,Rsx
