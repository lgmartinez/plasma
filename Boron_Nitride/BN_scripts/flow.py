import numpy
from parameters import wB,wNi,wCo,mB,mNi,mCo,mN,l2_target,l2_err,RSCD_N_B,RSCD_N_Ni,RSCD_N_Co
from initial_conditions import numx
from pressure_algorithm import get_pressure_GS
#--------------------------------------------------------
def get_density_bc(density,anode_bc_density):
    ''' updates densities in the anode region 
        by summing ablated density with local density
    
    Parameters:
    -----------
    density: total density in the domain
    anode_bc_density
    
    Returns:
    -------
    density with updated boundary values
        
    '''
    rhobc = density.copy()
    #anode region
    rhoan = anode_bc_density
    rhobc[0] = 2*rhoan - rhobc[1] #rhoan #
    #
    return rhobc
#-----------------------------------------------------------
def get_velocity_bc(velocity, ablation_velocity):
    ''' 
    Updates the velocity boundary conditions
    
    Parameters:
    -----------
    velocity: total (divergence-free) velocity
    ablation_velocity: velocity of sublimated particles
    
    Returns:
    ---------
    ubc: updated velocity boundary conditions
    '''
    ubc = velocity.copy()
    #anode region
    uxint0 = 2*(ablation_velocity) - 0.5*(ubc[2]+ablation_velocity)
    ubc[0] = 2*uxint0 - ablation_velocity 
    ubc[1] = ablation_velocity
    # cathode region
    ubc[-1] = 2*ubc[-2] - ubc[-3] #ubc[-3]# 
    #
    return ubc
#-----------------------------------------------------
def get_pressure_bc(pressure,anode_temperature,anode_density,anode_Rgas,\
                    cathode_temperature,cathode_density,cathode_Rgas):
    ''' 
    uses eqn of state to set pressure bc
    
    Parameters:
    -----------
    pressure: array 
    anode_temperature: first value in temperature array
    anode_density: first value in density array
    anode_Rgas: gas constant
    cathode_temperature 
    cathode_density
    cathode_Rgas

    Returns
    --------
    pbc: pressure with updated boundary conditions
    '''
    pbc = pressure.copy()
    pbc[0] = anode_temperature*anode_density*anode_Rgas
    #pbc[-1] = cathode_temperature*cathode_density*cathode_Rgas
    return pbc
#------------------------------------------------------
def get_interpolated_density(velocity,density, anode_density):
    '''
    Interpolates density values at the cell edges

    Returns:
    -------
    rhoxint: interpolated density array
    '''
    rhoan = anode_density
    rhox = density.copy()
    rhoxint = numpy.zeros_like(velocity)
    #ghost cell anode:
    rhoxint[0] = 2*rhox[0] - rhoan
    #anode-gap interface
    rhoxint[1] = rhoan
    #interior points:
    rhoxint[2:-1] = 0.5*(rhox[2:]+rhox[1:-1])
    #ghost cell cathode:
    rhoxint[-2] = 0.5*(rhox[-2]+rhox[-1]) 
    rhoxint[-1] = 2*rhoxint[-2] - rhoxint[-3]
    #
    return rhoxint
#----------------------------------------------------------
def get_momentum_bc(momentum,ablated_velocity,density_anode):
    ''' 
    Boundary Conditions for momentum
    
    Parameters:
    -----------
    momentum: momentum values
    ablated_velocity: anode bc velocity
    density_anode: anode bc density

    Returns:
    -------
    mombc: momentum with updated boundary conditions
    '''
    mombc = momentum.copy()
    moman = ablated_velocity*density_anode 
    # anode region
    mombc_int0 = 2*moman - 0.5*(moman+mombc[2])
    mombc[0] = 2*mombc_int0 - moman
    mombc[1] = moman
    # cathode/far-field region
    mombc[-1] = 2*mombc[-2] - mombc[-3]#mombc[-3]#
    #
    return mombc
#-----------------------------------------------------------
def get_convective_flux(velocity,density,delx,anode_density):
    ''' calculates convective fluxes for momentum equation
        uses central differencing scheme 
        
    Parameters:
    -----------
    velocity: helmholtz velocity
    delx: x-direction spacing
    
    Returns:
    -------
    fc: convective fluxes
    '''
    ux = velocity.copy()
    rhoxint = get_interpolated_density(velocity,density,anode_density)
    fc = numpy.zeros_like(velocity)
    #
    uxint = numpy.zeros_like(density)
    uxint[:] = 0.5*(ux[1:]+ux[:-1])
    uxint[-1] = 2*uxint[-2] - uxint[-3] #uxint[-2] #
    rhox = density.copy()
    #
    fc[1:-1] = (-0.5/delx)*((rhoxint[2:]*ux[2:]**2) -\
                            (rhoxint[:-2]*ux[:-2]**2) )
    # Alternative, use cell-center:
    fc[-2] = (-1./delx)*((rhox[-1]*uxint[-1]**2) -\
                            (rhox[-2]*uxint[-2]**2) )
    #
    return fc
#---------------------------------------------------------------
def get_viscous_flux(velocity,shear_viscosity,delx):
    ''' 
    calculates viscous fluxes for momentum equation
        
    Parameters:
    -----------
    velocity: helmholtz velocity
    delx: x-direction spacing
    shear_viscosity: mu
   
    Returns:
    --------
    fv: viscous fluxes
    '''
    ux = velocity.copy()
    fv = numpy.zeros_like(velocity)
    #
    mux = shear_viscosity.copy()
    muxb = numpy.zeros_like(velocity) #interpolated viscosity boundary values
    muxb[1:-1] = 0.5*(mux[1:]+mux[:-1])
    muxb[0] = muxb[1]
    muxb[-1] = muxb[-2]
    #
    uxint = numpy.zeros_like(shear_viscosity)
    uxint[:] = 0.5*(ux[1:]+ux[:-1])
    uxint[-1] = 2*uxint[-2] - uxint[-3] #uxint[-2]#
    #
    fv[1:-1] = (muxb[1:-1])*(1./(delx**2))*(ux[2:] - 2*ux[1:-1] + ux[:-2]) +\
                (1./delx)*(mux[1:]-mux[:-1])*\
                (0.5/delx)*(ux[2:]-ux[:-2])
    #
    fv[-2] = muxb[-2]*(2./(delx**2))*(uxint[-1]-2*ux[-2]+uxint[-2]) +\
            (1./delx)*(mux[-1]-mux[-2])*\
            (1./delx)*(uxint[-1]-uxint[-2])
    #
    return fv
#-------------------------------------------------------------------
def get_pressure_flux(velocity, pressure, delx):
    ''' 
    calculates pressure fluxes for momentum equation
    
    Parameters:
    ----------
    velocity: helmholtz velocity
    pressure: heavy particle pressure only, exclude electron pressure
    delx: x-direction spacing
    
    Returns:
    -------
    fp: pressure flux
    '''
    px = pressure.copy()
    fp = numpy.zeros_like(velocity)
    #
    fp[2:-2] = -(1./delx)*(px[2:-1]-px[1:-2])
    #
    return fp
#-------------------------------------------------------------------
def get_source_pressure(density_new, density_old, momentum, delx, delt):
    ''' 
    computes the pressure source terms
    
    Parameters:
    -----------
    velocity: helmholtz velocity
    pressure: heavy particle pressure
    delx: x-direction spacing
    
    Returns:
    --------
    fps: pressure source term
    '''
    rho_new = density_new.copy()
    rho_old = density_old.copy() 
    momx = momentum.copy()
    fps = numpy.zeros_like(density_old)
    #
    fps[1:-1] = (1./delx)*(momx[2:-1] - momx[1:-2]) +\
             (1./delt)*(rho_new[1:-1]-rho_old[1:-1])
    #
    return fps
#-------------------------------------------------------------------------
def get_flow(pressure,density,velocity,shear_viscosity,\
                 delx, delt, velocity_anode_bc, density_anode_bc):
    ''' Use Fractional Step Method to find velocity field and 
    pressure of heavy particles
    
    Parameters:
    -----------
    momentum:
    pressure:
    density:
    velocity: current velocity field
    shear_viscosity: mu
    delx: x-direction spacing
    delt: time step
    velocity_anode_bc
    density_anode_bc
    
    Returns:
    --------
    rhonew: new density values
    un: new velocity vaues
    ps3: new pressure values
    cflm: new momentum cfl number
    momfluxmax: maximum value of momentum flux
    '''
    # flux-related terms
    G1 = numpy.zeros_like(velocity)
    G2 = numpy.zeros_like(velocity)
    G3 = numpy.zeros_like(velocity)
    #Helmholtz momentum
    moms1 = numpy.zeros_like(velocity)
    moms2 = numpy.zeros_like(velocity)
    moms3 = numpy.zeros_like(velocity)
    # density arrays
    rhoold = density.copy()
    rho_bc_anode = density_anode_bc
    rhoxint = get_interpolated_density(velocity,density,density_anode_bc)
    # current velocity field
    un = velocity.copy()
    # current momentum
    momn = numpy.zeros_like(velocity)
    momn[:]=rhoxint[:]*un[:]
    # current shear viscosity field
    mun = shear_viscosity.copy()
    #anode bc velocity
    u_abl = velocity_anode_bc
    # current pressure
    pn = pressure.copy()
    anode_pressure = pn[0] # for pressure algorithm
    cathode_pressure = pn[-1] 
    #
    #Fractional Step Method, RK3 time advancement:
    #-------------------------------------------------------------------------
    #-------------------------------------------------------------------------
    #stage1
    #first thing is find the new density that would exist at this step
    rhonew = get_mass_conservation(rhoold,un,delx,(1./3.)*delt,\
                                   rho_bc_anode) #advance 1/3 dt
    #
    Fc1 = get_convective_flux(un,rhoold,delx,rho_bc_anode) #Convective Flux 
    Fv1 = get_viscous_flux(un,mun,delx) #Viscous FLux
    G1[:] = Fc1[:] + Fv1[:]
    #
    moms1[:] = momn[:] + (delt/3.)*(G1[:]) #predicted momentum
    moms1 = get_momentum_bc(moms1,u_abl,rho_bc_anode)
    #
    fs1 = get_source_pressure(rhonew,rhoold,moms1,delx,delt) #pressure source term 
    #
    ps1 = get_pressure_GS(pn,pn,delx,delt,l2_target,1./3.,fs1,numx,\
                          anode_pressure,cathode_pressure,l2_err)
    #
    Fp1 = get_pressure_flux(un, ps1, delx) #pressure flux
    #
    momn[:] = moms1[:] + (delt/3.)*(Fp1[:]) #rk1 momentum
    momn = get_momentum_bc(momn,u_abl,rho_bc_anode)
    #
    rhonew_int = get_interpolated_density(un,rhonew,rho_bc_anode)  
    un[:] = momn[:]/rhonew_int[:] #rk1 velocity
    un = get_velocity_bc(un, u_abl)
    #-------------------------------------------------------------------------
    #-------------------------------------------------------------------------
    #stage2
    rhoold = rhonew.copy()
    rhonew = get_mass_conservation(rhoold,un,delx,(5./12.)*delt,\
                                  rho_bc_anode) #advance to 9/12 dt
    #
    Fc2 = get_convective_flux(un,rhoold,delx,rho_bc_anode) #Convective Flux
    Fv2 = get_viscous_flux(un,mun,delx) #Viscous FLux
    G2[:] = - (5./9.)*(G1[:]) + (Fc2[:] + Fv2[:])
    #
    moms2[:] = momn[:] + (15.*delt/16.)*(G2[:]) 
    moms2 = get_momentum_bc(moms2,u_abl,rho_bc_anode)
    #
    fs2 = get_source_pressure(rhonew,rhoold,moms2,delx,delt) 
    #
    ps2 = get_pressure_GS(pn,pn,delx,delt,l2_target,5./12.,fs2,numx,\
                          anode_pressure,cathode_pressure,l2_err)
    #
    Fp2 = get_pressure_flux(un, ps2, delx) #pressure flux
    #
    momn[:] = moms2[:] + (5.*delt/12.)*(Fp2[:]) #rk2 momentum
    momn = get_momentum_bc(momn,u_abl,rho_bc_anode)
    #
    rhonew_int = get_interpolated_density(un,rhonew,rho_bc_anode) 
    un[:] = momn[:]/rhonew_int[:] #rk2 velocity
    un = get_velocity_bc(un, u_abl)
    #-------------------------------------------------------------------------
    #-------------------------------------------------------------------------
    #stage 3
    rhoold = rhonew.copy()
    rhonew = get_mass_conservation(rhoold,un,delx,(1./4.)*delt,\
                                  rho_bc_anode) #advance to dt
    #
    Fc3 = get_convective_flux(un,rhoold,delx,rho_bc_anode) #Convective Flux
    Fv3 = get_viscous_flux(un,mun,delx) #Viscous FLux
    G3[:] = - (153./128.)*(G2[:]) + (Fc3[:] + Fv3[:])
    #
    moms3[:] = momn[:] + (8.*delt/15.)*(G3[:]) #helmholtz velocity
    moms3 = get_momentum_bc(moms3,u_abl,rho_bc_anode)
    #
    fs3 = get_source_pressure(rhonew,rhoold,moms3,delx,delt) #source term for pressure
    #
    ps3 = get_pressure_GS(pn,pn,delx,delt,l2_target,1./4.,fs3,numx,\
                          anode_pressure,cathode_pressure,l2_err)
    #
    Fp3 = get_pressure_flux(un, ps3, delx) #pressure flux
    #
    momn[:] = moms3[:] + (1.*delt/4.)*(Fp3[:]) #rk1 velocity
    momn = get_momentum_bc(momn,u_abl,rho_bc_anode)
    #
    rhonew_int = get_interpolated_density(un,rhonew,rho_bc_anode) 
    un[:] = momn[:]/rhonew_int[:] #rk3 velocity
    un = get_velocity_bc(un, u_abl)
    #
    #momentum CFL:
    momfluxmax = numpy.max(rhonew_int[:]*un[:]**2)
    cflm = (delt/(delx))*momfluxmax
    #
    return rhonew, un, ps3, cflm,momfluxmax
#------------------------------------------------------------
def get_mass_conservation(density, velocity, delx, delt, anode_density):
    ''' solves mass conservation to obtain total density in the field
    
    Parameters:
    -----------
    density: current density distribution
    velocity: current velocity distribution
    delx: x-direction spacing
    delt: time step
    anode_density: input at anode-gap interface
    
    Returns:
    --------
    rho_new: new total density distribution
    
    '''
    # densities
    rhoan = anode_density
    rhox = density.copy()
    rho_new = numpy.zeros_like(density)
    #velocity
    ux = velocity.copy()
    ux_int = numpy.zeros_like(density)
    ux_int[:] = 0.5*(ux[1:]+ux[:-1])
    # mass conservation
    rho_new[1:-1] = rhox[1:-1] - (0.5*delt/delx)*\
                                    (rhox[2:]*ux_int[2:] - rhox[:-2]*ux_int[:-2])
    #boundary conditions
    rho_new[0] = 2*rhoan - rho_new[1] #rhoan #
    rho_new[-1] = 2*rho_new[-2] - rho_new[-3]#rho_new[-2] #
    #
    return rho_new
#-------------------------------------------------------------------
def get_mass_diffusion(pressure, temperature, density, velocity, concentration1,\
                      concentration2, concentration3, concentration4, delt, delx,\
                      evaporated_density, interface_anode_density):
    ''' Obtains individual species densities and concentrations by solving the
    mass diffusion equation
    
    Parameters:
    -----------
    pressure:
    temperature:
    density:
    velocity:
    concentration1:
    concentration2:
    concentration3
    concentration4
    delt:
    delx:
    evaporated_density:
    interface anode density:
    
    Returns:
    --------
    rho_s1: density of Boron
    rho_s2: density of Nickel
    rho_s3: density of Cobalt
    rho_g1: density of Nitrogen
    c1_out: concentration of Boron
    c2_out: concentration of Nickel
    c3_out: concentration of Cobalt
    c4_out: concentration of Nitrogen
    '''
    #temperature and pressure
    p = pressure.copy()
    Tx = temperature.copy()
    # current density
    rhox = density.copy()
    # densities of individual species
    rho_s1 = numpy.zeros_like(density) #B
    rho_s2 = numpy.zeros_like(density) #Ni
    rho_s3 = numpy.zeros_like(density) #Co
    rho_g1 = numpy.zeros_like(density) #N
    # velocity in domain
    ux = velocity.copy()
    uxc = numpy.zeros_like(density)
    uxc[:] = 0.5*(ux[1:]+ux[:-1])
    # concentrations
    c1 = concentration1.copy()
    c2 = concentration2.copy()
    c3 = concentration3.copy()
    cg1 = concentration4.copy()
    # new concentrations
    c1_out = numpy.zeros_like(concentration1)
    c2_out = numpy.zeros_like(concentration2)
    c3_out = numpy.zeros_like(concentration3)
    cg1_out = numpy.zeros_like(concentration4)
    #diffusion flux: D {d/dx rho*c}
    Dflux_s1 = numpy.zeros_like(density) 
    Dflux_s2 = numpy.zeros_like(density)
    Dflux_s3 = numpy.zeros_like(density)
    Dflux_g1 = numpy.zeros_like(density)
    # Rigid sphere collision diameters defined in parameters above
    # molar mass of base material and background gas
    ms1 = mB
    ms2 = mNi
    ms3 = mCo
    mg1 = mN
    #--------------------------------
    # diffusion coefficient
    D1 = numpy.zeros_like(density) #binary diffusion coefficient
    D1[1:-1] = (((2.63e-7)/(( p[1:-1]/101325.0)*(RSCD_N_B)**2))*\
                numpy.sqrt(((Tx[1:-1])**3)*(ms1+mg1)/(2*ms1*mg1)))*1e-4
    D2 = numpy.zeros_like(density) #binary diffusion coefficient
    D2[1:-1] = (((2.63e-7)/(( p[1:-1]/101325.0)*(RSCD_N_Ni)**2))*\
                numpy.sqrt(((Tx[1:-1])**3)*(ms2+mg1)/(2*ms2*mg1)))*1e-4
    D3 = numpy.zeros_like(density) #binary diffusion coefficient
    D3[1:-1] = (((2.63e-7)/(( p[1:-1]/101325.0)*(RSCD_N_Co)**2))*\
                numpy.sqrt(((Tx[1:-1])**3)*(ms3+mg1)/(2*ms3*mg1)))*1e-4
    #diffusion flux, species 1
    Dflux_s1[1:-1] = (D1[1:-1]/(2*delx))*(rhox[2:]*c1[2:] - rhox[:-2]*c1[:-2])
    Dflux_s1[0] = Dflux_s1[1]
    Dflux_s1[-1] = Dflux_s1[-2]
    #diffusion flux, species 2
    Dflux_s2[1:-1] = (D2[1:-1]/(2*delx))*(rhox[2:]*c2[2:] - rhox[:-2]*c2[:-2])
    Dflux_s2[0] = Dflux_s2[1]
    Dflux_s2[-1] = Dflux_s2[-2]
    #diffusion flux, species 3
    Dflux_s3[1:-1] = (D3[1:-1]/(2*delx))*(rhox[2:]*c3[2:] - rhox[:-2]*c3[:-2])
    Dflux_s3[0] = Dflux_s3[1]
    Dflux_s3[-1] = Dflux_s3[-2]
    # diffusion flux, gas
    Dflux_g1[1:-1] = (D1[1:-1]/(2*delx))*(rhox[2:]*cg1[2:] - rhox[:-2]*cg1[:-2])
    Dflux_g1[0] = Dflux_g1[1]
    Dflux_g1[-1] = Dflux_g1[-2]
    #---------------------------
    #species 1 density
    #anode conditions:
    rho_s1[0] = wB*evaporated_density + c1[1]*interface_anode_density #rhox[1]
    # inner points
    rho_s1[1:-1] = rhox[1:-1]*c1[1:-1] -\
                    (delt/(2*delx))*(rhox[2:]*uxc[2:]*c1[2:] -\
                                     rhox[:-2]*uxc[:-2]*c1[:-2]) +\
                    (delt/(2*delx))*(Dflux_s1[2:] - Dflux_s1[:-2])
    # exit
    rho_s1[-1] = 2*rho_s1[-2] - rho_s1[-3]#rho_s1[-2] #
    #----------------------------
    #species 2 density
    #anode conditions
    rho_s2[0] = wNi*evaporated_density + c2[1]*interface_anode_density #rhox[1]
    #
    # inner points
    rho_s2[1:-1] = rhox[1:-1]*c2[1:-1] -\
                    (delt/(2*delx))*(rhox[2:]*uxc[2:]*c2[2:] -\
                                     rhox[:-2]*uxc[:-2]*c2[:-2]) +\
                    (delt/(2*delx))*(Dflux_s2[2:] - Dflux_s2[:-2])
    # exit
    rho_s2[-1] = 2*rho_s2[-2] - rho_s2[-3]#rho_s2[-2] #
    #-----------------------------
    #----------------------------
    #species 3 density
    #anode conditions
    rho_s3[0] = wCo*evaporated_density + c3[1]*interface_anode_density #rhox[1]
    #
    # inner points
    rho_s3[1:-1] = rhox[1:-1]*c3[1:-1] -\
                    (delt/(2*delx))*(rhox[2:]*uxc[2:]*c3[2:] -\
                                     rhox[:-2]*uxc[:-2]*c3[:-2]) +\
                    (delt/(2*delx))*(Dflux_s3[2:] - Dflux_s3[:-2])
    # exit
    rho_s3[-1] = 2*rho_s3[-2] - rho_s3[-3]#rho_s3[-2]  #
    #-----------------------------
    #----------------------------
    #gas species density
    #anode conditions
    rho_g1[0] = cg1[1]*interface_anode_density #rhox[1]
    #
    # inner points
    rho_g1[1:-1] = rhox[1:-1]*cg1[1:-1] -\
                    (delt/(2*delx))*(rhox[2:]*uxc[2:]*cg1[2:] -\
                                     rhox[:-2]*uxc[:-2]*cg1[:-2]) +\
                    (delt/(2*delx))*(Dflux_g1[2:] - Dflux_g1[:-2])
    # exit
    rho_g1[-1] = 2*rho_g1[-2] - rho_g1[-3]#rho_g1[-2]  #
    #-----------------------------
    #concentration of species
    #
    rhototal = numpy.zeros_like(density)
    rhototal[:] = numpy.abs(rho_s1[:])+ numpy.abs(rho_s2[:]) +\
                   numpy.abs(rho_s3[:])+ numpy.abs(rho_g1[:])
    #
    c1_out[:] = numpy.abs(rho_s1[:])/rhototal[:]
    c2_out[:] = numpy.abs(rho_s2[:])/rhototal[:]
    c3_out[:] = numpy.abs(rho_s3[:])/rhototal[:]
    cg1_out[:] = numpy.abs(rho_g1[:])/rhototal[:]
    #
    return rho_s1, rho_s2,rho_s3, rho_g1, c1_out, c2_out, c3_out, cg1_out
