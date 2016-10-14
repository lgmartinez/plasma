import numpy
from parameters import molmass_evap,xB,xNi,xCo
'''
This script contains the following functions:

ablation.get_ablation_properties
'''
def get_ablation_properties(anode_temperature,local_pressure,\
                            local_temperature,local_Rgas):
    '''
    Parameters:
    -----------
    anode_temperature
    molar_fraction_s1: Boron
    molar_fraction_s2: Nickel
    molar_fraction_s3: Cobalt
    local_pressure: pressure at 1st node in gap
    local_temperature: temperature at 1st node in gap
    local_Rgas: gas constant at 1st node in gap

    Computes:
    ----------
    psat: total saturation pressure [Pa]
    abl_rate:  flux rate of evaporated mass [kg/(m2-s)]
    vth_evap: thermal velocity of evaporated mass [m/s]
    rho_evap: density of evaporated mass [kg/m3]
    rho_local: interface anode density (density of fluid in 1st node in gap)
    vel_anode: mass averaged velocity of sublimated particles (BC)
    rho_anode: total density incoming into domain (BC)

    Returns:
    --------
    rho_evap: sublimated density
    rho_local: density at interface
    rho_anode: rho_evap+rho_local
    vel_anode: mass-averaged bulk velocity of sublimated particles
    '''
    Ru = 8134. 
    #saturation pressure
    # Boron saturation pressure (using carbon curves for now)
    xs1 = xB
    A_s1 = 15.73
    B_s1 = 40030.0
    psat_s1 = (0.133*numpy.exp(2.3*(A_s1 - B_s1/anode_temperature)))
    # Nickel saturation pressure
    xs2 = xNi
    A_s2 = 12.75
    B_s2 = 20960.0
    psat_s2 = (0.133*numpy.exp(2.3*(A_s2 - B_s2/anode_temperature)))
    # Cobalt saturation pressure (use nickel for now)
    xs3 = xCo
    A_s3 = 12.75
    B_s3 = 20960.0
    psat_s3 = (0.133*numpy.exp(2.3*(A_s3 - B_s3/anode_temperature)))
    #total saturation pressure
    psat = xs1*psat_s1 + xs2*psat_s2 + xs3*psat_s3
    #
    #ablation rate
    abl_rate = psat*numpy.sqrt(molmass_evap/(2*numpy.pi*Ru*anode_temperature))
    #velocity of evaporated mass
    vth_evap = 0.5*numpy.sqrt((8*Ru*anode_temperature)/(molmass_evap*numpy.pi))
    #density of evaporated mass
    rho_evap = abl_rate/vth_evap
    #interface anode density
    rho_local = local_pressure/(local_temperature*local_Rgas)
    # mass-averaged bulk velocity, boundary condition
    vel_anode = abl_rate/(rho_evap+rho_local)
    #total anode density, used as boundary condition
    rho_anode = rho_evap+rho_local
    #done!
    return abl_rate,rho_evap,rho_local,rho_anode,vel_anode   
