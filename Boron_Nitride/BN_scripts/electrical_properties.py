#this script computes number densities, electron pressure
#this script also computes electrical conductivity
#Luis Martinez
#October 2, 2016
import scipy
from scipy import constants
import numpy
import math
from parameters import I_arc,Ips1,Ips2,Ips3,Ipg1,M_B,M_Ni,M_Co,M_N,Ran,Cme,Ce,Ckb,Ch

def saha_algorithm(temperature,density,concentration1,concentration2,concentration3,\
                   concentration4,electron_number_density):
    ''' 
    Solve Set of Saha equations to obtain number densities
    
    Parameters:
    -----------
    temperature
    density
    concentrations 1 through 4: Boron,Nickel,Cobalt,Nitrogen
    electron_number_density

    Returns:
    --------
    1. number densities of ions, neutrals, and electrons
    2. electron pressure
    '''
    #
    T_el = temperature.copy()
    ne_old = electron_number_density.copy()
    #
    rhox = density.copy()
    #
    c1 = concentration1.copy()
    c2 = concentration2.copy()
    c3 = concentration3.copy()
    cg1 = concentration4.copy()
    #
    E1 = Ips1
    E2 = Ips2
    E3 = Ips3
    Egas = Ipg1
    #
    M1 = M_B
    M2 = M_Ni
    M3 = M_Co
    Mg1 = M_N
    #
    nc1 = numpy.zeros_like(density)
    nc2 = numpy.zeros_like(density)
    nc3 = numpy.zeros_like(density)
    ncg1 = numpy.zeros_like(density)
    nc1[:] = rhox[:]*c1[:]/M1
    nc2[:] = rhox[:]*c2[:]/M2
    nc3[:] = rhox[:]*c3[:]/M3
    ncg1[:] = rhox[:]*cg1[:]/Mg1
    #
    SAHA_s1 = numpy.zeros_like(density)
    d1c1 = numpy.zeros_like(density)
    d1c1_old = numpy.zeros_like(density)
    as1 = numpy.zeros_like(density)
    bs1 = numpy.zeros_like(density)
    cs1 = numpy.zeros_like(density)
    d0c1 = numpy.zeros_like(density)
    de_c1 = numpy.zeros_like(density)
    s1_ions = numpy.zeros_like(density)
    err1 = numpy.zeros_like(density)
    #
    SAHA_s2 = numpy.zeros_like(density)
    d1c2 = numpy.zeros_like(density)
    d1c2_old = numpy.zeros_like(density)
    as2 = numpy.zeros_like(density)
    bs2 = numpy.zeros_like(density)
    cs2 = numpy.zeros_like(density)
    d0c2 = numpy.zeros_like(density)
    de_c2 = numpy.zeros_like(density)
    s2_ions = numpy.zeros_like(density)
    err2 = numpy.zeros_like(density)
    #
    SAHA_s3 = numpy.zeros_like(density)
    d1c3 = numpy.zeros_like(density)
    d1c3_old = numpy.zeros_like(density)
    as3 = numpy.zeros_like(density)
    bs3 = numpy.zeros_like(density)
    cs3 = numpy.zeros_like(density)
    d0c3 = numpy.zeros_like(density)
    de_c3 = numpy.zeros_like(density)
    s3_ions = numpy.zeros_like(density)
    err3 = numpy.zeros_like(density)
    #
    SAHA_gas = numpy.zeros_like(density)
    d1g = numpy.zeros_like(density)
    d1g_old = numpy.zeros_like(density)
    ag = numpy.zeros_like(density)
    bg = numpy.zeros_like(density)
    cg = numpy.zeros_like(density)
    d0g = numpy.zeros_like(density)
    de_g = numpy.zeros_like(density)
    gas_ions = numpy.zeros_like(density)
    errg = numpy.zeros_like(density)
    ##
    s1_neutrals = numpy.zeros_like(density)
    s2_neutrals = numpy.zeros_like(density)
    s3_neutrals = numpy.zeros_like(density)
    gas_neutrals = numpy.zeros_like(density)
    electrons = numpy.zeros_like(density)
    P_electrons = numpy.zeros_like(density) #electron pressure
    #
    ions_total = numpy.zeros_like(density)
    neutrals_total = numpy.zeros_like(density)
    #
    d1c1[:] = 0.0
    d1c2[:] = 0.0
    d1c3[:] = 0.0
    d1g[:] = 0.0
    #
    err1[:] = 10.
    err2[:] = 10.
    err3[:] = 10.
    errg[:] = 10.
    #
    d0c1[:] = 1.0
    d0c2[:] = 1.0
    d0c3[:] = 1.0
    d0g[:] = 1.0
    #
    iterations = 0
    itermax = 20000
    itermax2 = itermax-10
    tolerance = 1e-7
    #
    domain_size = density.size
    #
    for i in range(domain_size):
        #
        while ((err1[i] > tolerance) and (err2[i]>tolerance) and\
               (err3[i]>tolerance) and (errg[i]>tolerance)):
            #----------------------------------------------------------------------
            #species 1
            d1c1_old[i] = d1c1[i]
            SAHA_s1[i] =  ( ( (2*numpy.pi*Cme*Ckb*T_el[i]) / (Ch**2) )**(3/2) ) *\
                     math.exp((-E1)/(constants.k*T_el[i]))
            as1[i] = nc1[i]
            bs1[i] = nc2[i]*d1c2[i] + nc3[i]*d1c3[i] + ncg1[i]*d1g[i] + SAHA_s1[i]
            cs1[i] = - SAHA_s1[i]
            d1c1[i] = (-bs1[i] + numpy.sqrt(bs1[i]**2 - 4*as1[i]*cs1[i]))/(2*as1[i])
            #
            d0c1[i] = 1. - d1c1[i]
            de_c1[i] = d1c1[i]
            #
            s1_ions[i] = d1c1[i]*rhox[i]*c1[i]/M1
            #-----------------------------------------------------------------
            #species 2
            d1c2_old[i] = d1c2[i]
            SAHA_s2[i] =  ( ( (2*numpy.pi*Cme*Ckb*T_el[i]) / (Ch**2) )**(3/2) ) *\
                     math.exp((-E2)/(constants.k*T_el[i]))
            as2[i] = nc2[i]
            bs2[i] = nc1[i]*d1c1[i] + nc3[i]*d1c3[i] + ncg1[i]*d1g[i] + SAHA_s2[i]
            cs2[i] = - SAHA_s2[i]
            d1c2[i] = (-bs2[i] + numpy.sqrt(bs2[i]**2 - 4*as2[i]*cs2[i]))/(2*as2[i])
            #
            d0c2[i] = 1. - d1c2[i]
            de_c2[i] = d1c2[i]
            #
            s2_ions[i] = d1c2[i]*rhox[i]*c2[i]/M2
            #-----------------------------------------------------------------
            #species 3
            d1c3_old[i] = d1c3[i]
            SAHA_s3[i] =  ( ( (2*numpy.pi*Cme*Ckb*T_el[i]) / (Ch**2) )**(3/2) ) *\
                     math.exp((-E3)/(constants.k*T_el[i]))
            as3[i] = nc3[i]
            bs3[i] = nc1[i]*d1c1[i] + nc2[i]*d1c2[i] + ncg1[i]*d1g[i] + SAHA_s3[i]
            cs3[i] = - SAHA_s3[i]
            d1c3[i] = (-bs3[i] + numpy.sqrt(bs3[i]**2 - 4*as3[i]*cs3[i]))/(2*as3[i])
            #
            d0c3[i] = 1. - d1c3[i]
            de_c3[i] = d1c3[i]
            #
            s3_ions[i] = d1c3[i]*rhox[i]*c3[i]/M3
            #---------------------------------------------------------------------
            #gas
            d1g_old[i] = d1g[i]
            SAHA_gas[i] = ( ( (2*numpy.pi*Cme*Ckb*T_el[i]) / (Ch**2) )**(3/2) ) *\
                     math.exp((-Egas)/(constants.k*T_el[i]))
            ag[i] = ncg1[i]
            bg[i] = nc1[i]*d1c1[i] + nc2[i]*d1c2[i] + nc3[i]*d1c3[i] + SAHA_gas[i]
            cg[i] = - SAHA_gas[i]
            d1g[i]  = (-bg[i] + numpy.sqrt(bg[i]**2 - 4*ag[i]*cg[i]))/(2*ag[i])
            #
            d0g[i] = 1. - d1g[i]
            de_g[i] = d1g[i]
            #
            gas_ions[i] = d1g[i]*rhox[i]*cg1[i]/Mg1
            #-------------------------------------------------------------------
            # electrons
            electrons[i] = s1_ions[i] + s2_ions[i] +\
                           s3_ions[i] + gas_ions[i] #number density, enforce LTE
            P_electrons[i] = numpy.abs(electrons[i]-ne_old[i])*Ckb*T_el[i]
            #--------------------------------------------------------------------
            #neutrals
            s1_neutrals[i] = d0c1[i]*rhox[i]*c1[i]/M1
            s2_neutrals[i] = d0c2[i]*rhox[i]*c2[i]/M2
            s3_neutrals[i] = d0c3[i]*rhox[i]*c3[i]/M3
            gas_neutrals[i] = d0g[i]*rhox[i]*cg1[i]/Mg1
            #
            #--------------------------------------------------------------------
            #totals
            ions_total[i] = gas_ions[i] + s1_ions[i] + s2_ions[i] + s3_ions[i] 
            neutrals_total[i] = s1_neutrals[i] + s2_neutrals[i] +\
                                s3_neutrals[i] + gas_neutrals[i]
            #----------------------------------------------------------------------
            #bound for minimum number of electrons:
            if (electrons[i]<ne_old[i]):
                electrons[i] = ne_old[i]
            if (ions_total[i]<ne_old[i]):
                ions_total[i] = ne_old[i]#electrons[i]
            #
            #errors
            err1[i] = numpy.abs(d1c1[i] - d1c1_old[i])
            err2[i] = numpy.abs(d1c2[i] - d1c2_old[i])
            err3[i] = numpy.abs(d1c3[i] - d1c3_old[i])
            errg[i] = numpy.abs(d1g[i] - d1g_old[i])
            #
            iterations += 1
        
            if (iterations>itermax):
                print('max iterations reached')
                break
            #
    return electrons,ions_total, neutrals_total,P_electrons
#------------------------------------------------------------------------
def get_electrical_conductivity(temperature,electron_number_density,\
                                neutral_number_density,time):
    ''' 
    Obtain electrical conductivity based on the chapmann-enskog equation
    
    Parameters:
    ------
    temperature: temperature obtained from energy equation
    electron_number_density
    neutral_number_density
    time
    
    Return:
    -------
    e_cond: electrical conductivitiy
    '''
    domain_size = numpy.size(temperature)
    Qm = numpy.zeros((domain_size),dtype=float)
    #
    na = neutral_number_density.copy()
    ne = electron_number_density.copy()
    #
    Tx = temperature.copy()
    #
    #momentum transfer cross section:
    r1 = 5.25
    r2 = 7.0
    Qm[1:-1] = ((r1 + (r2-r1)*(Tx[1:-1]-300.)/(11604.*(1.74-0.01)))*1e-16)*1e-4
    Qm[0] = Qm[1]
    Qm[-1] = Qm[-2]
    #Qm[:] = 40e-20
    #
    #electron-neutral collision
    v_e_a = numpy.zeros((domain_size),dtype=float)
    v_e_a[1:-1] = (4./3.)*Qm[1:-1]*(na[1:-1])*\
                    numpy.sqrt((8.*Ckb*Tx[1:-1]) /(numpy.pi*Cme))
    v_e_a[0] = v_e_a[1]
    v_e_a[-1] = v_e_a[-2]
    #electron-ion collision
    v_e_i = numpy.zeros((domain_size),dtype=float)
    lnA = numpy.zeros((domain_size),dtype=float)
    ke = numpy.zeros((domain_size),dtype=float)
    #
    ke[1:-1] = numpy.sqrt((4.*numpy.pi*ne[1:-1]*Ce**2)/(Ckb*Tx[1:-1]))
    gam = numpy.exp(0.577)
    lnA[1:-1] = numpy.log(4*Ckb*Tx[1:-1]/(ke[1:-1]*(Ce**2)*(gam**2))) -\
                2*numpy.log(numpy.sqrt(2.))
    #
    v_e_i[1:-1] = lnA[1:-1]*(4./3.)*(numpy.sqrt(2*numpy.pi))*ne[1:-1]*\
                    numpy.sqrt((Ckb*Tx[1:-1])/Cme)*\
                    (Ce**2/(Ckb*Tx[1:-1]))**2
    v_e_i[0] = v_e_i[1]
    v_e_i[-1] = v_e_i[-2]
    #electrical conductivity
    #e_cond =  numpy.zeros((domain_size),dtype=float)
    e_cond_real =  numpy.zeros((domain_size),dtype=float)
    econd_max = 2e4
    econd_min = I_arc
    #
    e_cond_real[1:-1] = ((Ce**2)/Cme)*(ne[1:-1]/(v_e_i[1:-1]+v_e_a[1:-1]))
    e_cond_real[0] = e_cond_real[1]
    e_cond_real[-1] = e_cond_real[-2]
    #
    for i in range(domain_size):
        if e_cond_real[i]>econd_max:
            e_cond_real[i] = econd_max
        if e_cond_real[i]<econd_min:
            e_cond_real[i] = econd_min
    #
    #if too many fluctuations, can set the entire domain to the average:
    econdavg = numpy.zeros((domain_size),dtype=float)
    econdavg = numpy.average(e_cond_real[:])
    e_cond_real[:] = econdavg
    #
    # e_cond[:] = econd_min
    #
    return e_cond_real
#-----------------------------------------------------------------------------
from parameters import Ran
def get_electric_potential(current_density,electrical_conductivity,\
                           electric_potential, center_nodes, delx):
    ''' Defines the electric potential field.
    
    This function obtains the electric potential based on two solutions to the 
    ordinary differential equation.
    
    Solution 1: B^2 - 4AC = 0
    Solution 2: B^2 - 4AC > 0
    
    Parameters
    ----------
    current_density: current density in the domain, constant current is assumed
    electrical_conductivity: 
    electric_potential: electric potential is solved from d/dx{-econd d/dx(phi)} =0
    center_nodes: cell-center nodes are used to solve electric potential
    delx: x-direction spacing
    
    Returns:
    --------
    electric potential
    '''
    jx = current_density.copy()
    #
    xc_nodes = center_nodes.copy()
    x_an = xc_nodes[0] - 0.5*delx
    x_cath = xc_nodes[-1] + 0.5*delx
    #
    array_size = electrical_conductivity.size
    eratio = numpy.zeros_like(center_nodes)
    #
    A = numpy.zeros_like(electrical_conductivity)
    B = numpy.zeros_like(electrical_conductivity)
    #
    C1 = numpy.zeros_like(electrical_conductivity)
    C2 = numpy.zeros_like(electrical_conductivity)
    #
    econdx = electrical_conductivity.copy()
    phix = electric_potential.copy()
    #
    jx_anode = I_arc/(numpy.pi*Ran**2) 
    econdx_anode =2e4 #electrical conductivity max
    anode_potential = jx_anode*delx/econdx_anode
    #
    cathode_potential =  -1.
    #
    A[1:-1] =  - econdx[1:-1]
    B[1:-1] = - (0.5/delx)*(econdx[2:] - econdx[:-2])
    #
    for i in range(1,array_size-1):
        if (B[i]**2 == 0.0 ):
            # general solution constants
            C2[i] = (cathode_potential/(numpy.exp((-B[i]/A[i])*x_cath)) -\
                     anode_potential)/(x_cath-x_an)
            C1[i] = (anode_potential/(numpy.exp((-B[i]/A[i])*x_an))) -\
                     anode_potential
            # potential solution
            phix[i] = C1[i]*numpy.exp((-B[i]/A[i])*xc_nodes[i-1]) +\
                      C2[i]*xc_nodes[i-1]*numpy.exp((-B[i]/A[i])*xc_nodes[i-1])
        if (B[i]**2 > 0.0):
            #print('using {B**2 > 0} solution for potential')
            eratio[:] = numpy.exp((-B[i]/A[i])*x_cath)/numpy.exp((-B[i]/A[i])*x_an)
            C2[i] = (cathode_potential-anode_potential*eratio[i-1])/(1.-eratio[i-1])
            C1[i] = (anode_potential-C2[i])/numpy.exp((-B[i]/A[i])*x_an)
            #potential solution
            phix[i] = C1[i]*numpy.exp((-B[i]/A[i])*xc_nodes[i-1]) + C2[i]
    #
    phix[0] = jx_anode*delx/econdx_anode
    phix[-1] = cathode_potential
    #
    return phix
