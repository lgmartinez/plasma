import numpy
from parameters import CAG,Ce,Ckb,anode_workfunction,cathode_workfunction,Cme,dH_an,cp_an,k_an,rho_an_max,k_cath,rho_cath_max,cp_cath,Uiz_cath
#
def get_energy(ablation_rate,density,velocity,pressure,reference_pressure,temperature,\
               enthalpy,specific_heat,thermal_conductivity,viscosity,\
               current_density,electrical_conductivity,electron_density,delx,delt):
    ''' 
    Defines energy based on computed values, 
    includes Joule Heating, Heat Flux, DP/Dt, d(puh)/dx
    
    Note that this actually solves for (rho*h), so 
                                         h = (rho*h)\rho
    
    Then we can solve for temperature,
                                         T^{n+1} = dh/cp + T_{ref}
    
    Parameters:
    ----------
    ablation_rate
    density
    velocity
    pressure
    reference_pressure
    temperature
    enthalpy
    specific_heat
    thermal_conductivity
    viscosity
    current_density
    electrical_conductivity
    electron_density
    delx
    delt

    Returns:
    --------
    Temperature at next time step
    enthalpy at next time step
    cfl no.
    max flux value
    anode heat flux total
    '''
    #density
    rhox = density.copy() #n
    #velocity
    ux = velocity.copy() #n
    uxint = numpy.zeros_like(enthalpy)
    uxint[:] = 0.5*(ux[1:]+ux[:-1])
    #pressure
    P_ref = reference_pressure.copy()#pressure.copy() # #n-1
    Px = pressure.copy() #n
    #temperature
    Tx = temperature.copy() #n
    #enthalpy
    hx = enthalpy.copy() # n
    #thermophysical properties
    cpx = specific_heat.copy()
    kx = thermal_conductivity.copy()
    mux = viscosity.copy()
    #electrical properties
    jx = current_density.copy()
    econdx = electrical_conductivity.copy()
    #electron density
    nex = electron_density.copy()
    #
    # create flux arrays for energy equation
    Vflux = numpy.zeros_like(enthalpy) #viscous heat flux
    Hflux = numpy.zeros_like(enthalpy) #conductive heat flux
    Jflux = numpy.zeros_like(enthalpy) #joule heating
    Pflux = numpy.zeros_like(enthalpy) # pressure flux
    Mflux  = numpy.zeros_like(enthalpy) #momentum flux
    rhoh_old = numpy.zeros_like(enthalpy) # rho*h @ n-1
    #property arrays
    rhoh_new = numpy.zeros_like(enthalpy) #density*enthalpy n+1
    h_new = numpy.zeros_like(enthalpy) #enthalpy n+1
    T_new = numpy.zeros_like(enthalpy) #temperature n+1
    Temp = numpy.zeros_like(enthalpy) # delta Temp
    Temp[:] = hx[:]/cpx[:] #Tx[:]#
    #
    # joule heating flux
    Jflux[1:-1] = delt*((jx[1:-1]**2)/econdx[1:-1])+\
                    delt*((5./2.)*Ckb/Ce)*(jx[1:-1])*\
                    (0.5/delx)*(Temp[2:] - Temp[:-2])
    #viscous heat flux
    vel_mu = numpy.zeros_like(enthalpy)
    vel_mu[1:-1] = uxint[1:-1]*mux[1:-1]
    vel_mu[0] = vel_mu[1]
    vel_mu[-1] = vel_mu[-2]
    Vflux[1:-1] = -(0.5*delt/delx)*(vel_mu[2:]*uxint[2:]-vel_mu[:-2]*uxint[:-2]) -\
                -(vel_mu[1:-1])*(0.5*delt/delx)*(uxint[2:]-2*uxint[1:-1]-uxint[:-2])
    # conductive heat flux
    Hflux[1:-1] = delt*(kx[1:-1]/(delx**2))*(Temp[2:] - 2*Temp[1:-1] + Temp[:-2])+\
                  delt*(0.5/delx)*((kx[2:]-kx[:-2])*(Temp[2:]-Temp[:-2]))
    # pressure flux
#     Pflux[1:-1] = (Px[1:-1]-P_ref[1:-1]) +\
#                     (0.5*delt/delx)*(Px[2:]-Px[:-2])*uxint[1:-1]
    Pflux[1:-1] = -(0.5*delt/delx)*(Px[2:]*uxint[2:]-Px[:-2]*uxint[:-2])
    # convective heat flux, or momentum flux
    Mflux[1:-1] = - (0.5*delt/delx)*\
                    (rhox[2:]*uxint[2:]*hx[2:] - rhox[:-2]*uxint[:-2]*hx[:-2])
    # Energy at n:
    rhoh_old[1:-1] = rhox[1:-1]*hx[1:-1]
    # Energy at n+1
    rhoh_new[1:-1] = rhoh_old[1:-1] + Mflux[1:-1] + Pflux[1:-1] +\
                        Hflux[1:-1] + Jflux[1:-1]
    # new (n+1) enthalpy values at interior points
    #note that here we approximate by dividing with rho @ n bc we dont have @ n+1
    h_new[1:-1] = rhoh_new[1:-1]/(rhox[1:-1])
    #------------------------------------------------
    # new temperature values
    #interior points
    T_new[1:-1] = ((h_new[1:-1]-hx[1:-1])/cpx[1:-1]) + Tx[1:-1]
    #T_new[1:-1] = ((h_new[1:-1])/cpx[1:-1])
    T_new[0] = T_new[1]
    T_new[-1] = T_new[-2]
    #temperature boundary conditions
    T_new,Tadd_a,Tadd_c = get_temperature_bc(T_new,nex[1],jx[1],ablation_rate,\
                                             jx[-2],nex[-2],delx,delt)
    #enthalpy boundary values
    h_new[0] = cpx[0]*T_new[0]
    h_new[-1] = h_new[-2] #cpx[-1]*T_new[-1]
    # maximum flux and cfl number
    efluxmax = numpy.max(rhoh_new[:]*numpy.abs(uxint[:]))
    cfle = (delt/(delx))*efluxmax
    #
    return T_new,h_new,cfle,efluxmax,Tadd_a,Tadd_c
#-------------------------------------------------------------------------------------
def get_temperature_bc(temperature,electron_density_anode,\
                       current_density_anode,ablation_rate,\
                       electron_density_cathode,current_density_cathode,delx,delt):
    ''' 
    Obtains boundary temperature at anode from heat fluxes
    
    parameters:
    ----------
    temperature: temperature at current time step
    electron_density_anode: use index=1
    current_density_anode: anode current density
    ablation_rate
    electron_density_cathode
    current_density_cathode
    delx
    delt
    
    Returns:
    --------
    temperature with boundary conditions
    '''
    Tbc = temperature.copy()
    Tan_ref = Tbc[0]
    #
    ne_an = electron_density_anode 
    j_an = current_density_anode 
    #
    abl_flux = ablation_rate
    #--------------------ANODE REGION-------------------------------------------#
    #work function of anode (eV)
    phi_an_eV = anode_workfunction 
    #temperature of oncoming electrons in kevlin
    Te_an = Tbc[1] 
    #temperature of oncoming electrons in eV
    Tev_an = Te_an/(11604.52500617)
    #electron flux [1/(m^2 sec)]=vth*ne
    jth = Ce*ne_an*numpy.sqrt((Ckb*Te_an)/(2*numpy.pi*Cme)) #j_an#
    #anode sheath drop, eV
    Ua_eV = - Tev_an*numpy.log(jth/j_an) 
    #anode heat flux, W/m^2 (oer unit length)
    qa = j_an*numpy.abs(2*Tev_an + phi_an_eV + Ua_eV)  
    # Heat of vaporization of anode material flux
    Levap_flux = dH_an*abl_flux
    # enthalpy of evaporated mass flux
    #hevap_flux = cp_an*(Tbc[0])*uan*rhoan 
    # volumetric heat flux, anode (per unit length) 
    heat_flux_anode = qa - Levap_flux 
    # define anode temperature
    Tan_add = (delt/(rho_an_max*cp_an))*\
              ((0.5*k_an/(delx**2))*(Tbc[2]-Tbc[0]) + heat_flux_anode)
    Tbc[0] = Tan_ref + Tan_add
    #
    #if (heat_flux_anode <= 0.):
       # Tbc[0] = Tbc[1]
    #--------------------Cathode REGION-------------------------------------------#
    Tbc[-1] = 2*Tbc[-2] - Tbc[-3] #Tbc[-2] #
    #Tcath_ref = Tbc[-1]
    #ne_cath = electron_density_cathode
    j_cath = current_density_cathode
    phi_cath_eV = cathode_workfunction
    #Te_cath = Tbc[-2]
    #Tev_cath =  Te_cath/(11604.52500617)
    #jth_cath = Ce*ne_cath*numpy.sqrt((Ckb*Te_cath)/(2*numpy.pi*Cme))
    #Uc_eV = -Tev_cath*numpy.log(jth_cath/j_cath)
    #qc = j_cath*(2*Tev_cath + phi_cath_eV  + Uc_eV )
    phi_cath_J = 1.60218e-19*cathode_workfunction
    j_em = CAG*(Tbc[-2]**2)*numpy.exp(phi_cath_J/(Ckb*Tbc[-2]))
    # include emission current only after this condition is satisfied:
    if (Tbc[-2]<3800):
        qc = j_cath*(phi_cath_eV + Uiz_cath)
    if (Tbc[-2]>= 3800):
        qc = j_cath*(phi_cath_eV + Uiz_cath) - j_em*phi_cath_eV
    qc = j_cath*(phi_cath_eV + Uiz_cath) - j_em*phi_cath_eV
    cathode_heat_flux = qc # W/m2 per unit length
    Tcath_add = (delt/(rho_cath_max*cp_cath))*\
                ((0.5*k_cath/(delx**2))*(Tbc[-3]-Tbc[-1]) + cathode_heat_flux)
#    Tbc[-1] = Tcath_ref + Tcath_add
    #
    return Tbc,Tan_add,Tcath_add
