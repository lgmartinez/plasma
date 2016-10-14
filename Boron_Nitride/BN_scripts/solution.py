import numpy
from parameters import mN,mB,mNi,mCo,k_an,cp_an,l2_target
from thermophysical_properties import get_thermophysical_properties
from time_step import get_adaptive_time_step
from ablation import get_ablation_properties
from initial_conditions import numx,xc,xb
from flow import get_velocity_bc, get_pressure_bc, get_density_bc,\
                 get_flow, get_mass_diffusion
from compressibility import get_compressibility
from electrical_properties import saha_algorithm,get_electrical_conductivity
from energy import get_energy
from plot_functions import plot, plot2
from matplotlib import pyplot
import matplotlib
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
#
def get_solution(enthalpy,temperature,reference_temperature,pressure,\
                 reference_pressure,density,reference_density,velocity,\
                 concentration_s1,concentration_s2,concentration_s3,\
                 concentration_g1,electric_potential,current_density,\
                 electrical_conductivity,number_density,\
                 delx,delt,cfl_energy,cfl_momentum,energyfluxmax,\
                 momentumfluxmax,reporting_interval1,reporting_interval2):
    ''' 
    this function combines the entire algorithm and provides a
    time-marched solution

    Parameters:
    -----------
    enhtalpy
    temperature
    reference_temperature
    pressure
    reference_pressure
    density
    reference_density
    velocity
    concentrations 1 through 4: B,Ni.Co,N
    electric_potential
    current_density
    electrical_conductivity
    number_density (e)
    delx,delt
    cfl_energy,cfl_momentum (if no initial value: use 125,0.75 respectively)
    energyfluxmax,momentumfluxmax (if no inital value: use 1e5,1.0 respectively)
    reporting_interval1,reporting_interval2

    Returns:
    -----------
    for now, output 2 scalar variables:
    number of iterations
    final time
    '''
    hnref = enthalpy.copy()
    #used for stability condition
    delt_ref = delt
    # cfl numbers for stability condition
    cflen = cfl_energy
    cflmn = cfl_momentum
    cflen_ref = cfl_energy
    cflmn_ref = cfl_momentum
    #fluxes for stability condition
    efmax = energyfluxmax
    mfmax = momentumfluxmax
    #temperature
    Tnref = temperature.copy()
    Tnref_old = reference_temperature.copy()
    #pressure
    pxn = pressure.copy() 
    pxn_ref = reference_pressure.copy()
    #density
    rhoxn = density.copy()
    rhoxn_ref = reference_density.copy()
    #velocity
    uxn = velocity.copy()
    #concentrations
    cB  = concentration_s1.copy()
    cNi = concentration_s2.copy()
    cCo = concentration_s3.copy()
    cN  = concentration_g1.copy() 
    # electric potential
    phin = electric_potential.copy()
    #current density
    jn = current_density.copy()
    # electrical conductivity
    econdxn = electrical_conductivity.copy()
    #number density of electrons
    nexn = number_density.copy()
    nexn_ref = number_density.copy() #reference_number_density
    # current time
    time = 0 + delt_ref
    #final time in seconds
    t_terminal = 0.0005
    #iteration count
    iterations = 0
    #-------------------
    #file = open('output_summary.txt','a')
    #
    while (time <t_terminal):
        iterations += 1
        delt = get_adaptive_time_step(cflen,cflmn,efmax,mfmax,\
                                      delx,delt_ref,numpy.max(uxn))
        time = time + delt
        mun,ktn,cpn,Rsn = get_thermophysical_properties(Tnref,cN,cB,cNi,cCo)
        ablrate,rhoevap,rholocal,rhoanbc,uanbc = get_ablation_properties(Tnref[0],pxn[1],\
                                                                 Tnref[1],Rsn[1])
        uxn = get_velocity_bc(uxn,uanbc)
        rhoxn = get_density_bc(rhoxn,rhoanbc)
        pxn = get_pressure_bc(pxn,Tnref[0],rhoxn[0],Rsn[0],\
                              Tnref[-1],rhoxn[-1],Rsn[-1]) 
        rhoxn,uxn,pxn,cflmn,mflux = get_flow(pxn,rhoxn,uxn,mun,\
                                             delx,delt,uanbc,rhoanbc)
        compress = get_compressibility(pxn,pxn_ref,rhoxn,rhoxn_ref)
        rhoB,rhoNi,rhoCo,rhoN,\
        cB,cNi,cCo,cN = get_mass_diffusion(pxn,Tnref,rhoxn,uxn,cB,cNi,cCo,cN,\
                                           delt,delx,rhoevap,rholocal)
        nexn,nixn,noxn,pexn = saha_algorithm(Tnref,rhoxn,cB,cNi,cCo,cN,nexn_ref)
        econdxn = get_electrical_conductivity(Tnref,nexn,noxn,time)
        Tn,hn,cflen,efmax,Tpa,Tpc = get_energy(ablrate,rhoxn,uxn,pxn,pxn_ref,Tnref,hnref,cpn,\
                                           ktn,mun,jn,econdxn,nexn,delx,delt)
        pxn[:] = rhoxn[:]*Rsn[:]*Tn[:]
	#reference values
        Tnref_old = Tnref.copy()
        Tnref = Tn.copy()
        hnref = hn.copy()
        rhoxn_ref = rhoxn.copy()
        pxn_ref = pxn.copy()
        nexn_ref = nexn.copy()
        delt_ref = delt
        time_ref = time
        cflen_ref = cflen
        cflmn_ref = cflmn
	#print summary
        if iterations in reporting_interval1:
            file = open('output_summary.txt','a')
            file.write('---------------------------------------\n' )
            file.write('current iteration is: %.3g\n' %iterations)
            file.write('evaporated density is: %.3g\n' %rhoevap)
            file.write('ablation rate is: %.3g\n' %ablrate)
            file.write('max compressibility is: %.3g \n' %numpy.max(compress))
            file.write('time step value is: %.3g\n' %delt)
            file.write('current time in seconds is: %.3g\n' %time)
            file.write('uxn[0] velocity in [m/s] is: %.3g\n' %uxn[0])
            file.write('uxn[1] velocity in [m/s] is: %.3g\n' %uxn[1])
            file.write('uxn[2] velocity in [m/s] is: %.3g\n' %uxn[2])
            file.write('uxn[3] velocity in [m/s] is: %.3g\n' %uxn[3])
            file.write('uxn[-1] velocity in [m/s] is: %.3g\n' %uxn[-1])
            file.write('uxn[-2] velocity in [m/s] is: %.3g\n' %uxn[-2])
            file.write('uxn[-3] velocity in [m/s] is: %.3g\n' %uxn[-3])
            file.write('max velocity is:  %.3g\n' %numpy.max(uxn))
            file.write('min velocity is:  %.3g\n' %numpy.min(uxn))
            file.write('anode density in [kg/m3] is: %.3g\n' %rhoxn[0])
            file.write('rhoxn[1] density in [kg/m3] is: %.3g\n' %rhoxn[1])
            file.write('rhoxn[2] density in [kg/m3] is: %.3g\n' %rhoxn[2])
            file.write('rhoxn[3] density in [kg/m3] is: %.3g\n' %rhoxn[3])
            file.write('rhoxn[-1] density in [kg/m3] is: %.3g\n' %rhoxn[-1])
            file.write('rhoxn[-2] density in [kg/m3] is: %.3g\n' %rhoxn[-2])
            file.write('rhoxn[-3] density in [kg/m3] is: %.3g\n' %rhoxn[-3])
            file.write('max density is:  %.3g\n' %numpy.max(rhoxn))
            file.write('min density is:  %.3g\n' %numpy.min(rhoxn))
            file.write('Anode Temperature [K] is: %.3g\n' %Tn[0])
            file.write('T[1] temperature in [K] is: %.3g\n' %Tn[1])
            file.write('T[2] temperature in [K] is: %.3g\n' %Tn[2])
            file.write('T[3] temperature in [K] is: %.3g\n' %Tn[3])
            file.write('T[-1] temperature in [K] is: %.3g\n' %Tn[-1])
            file.write('T[-2] temperature in [K] is: %.3g\n' %Tn[-2])
            file.write('T[-3] temperature in [K] is: %.3g\n' %Tn[-3])
            file.write('max temperature is:  %.3g\n' %numpy.max(Tn))
            file.write('min temperature is:  %.3g\n' %numpy.min(Tn))
            file.write('Anode Pressure [Pa] is: %.3g\n' %pxn[0])
            file.write('p[1] pressure in [Pa] is: %.3g\n' %pxn[1])
            file.write('p[2] pressure in [Pa] is: %.3g\n' %pxn[2])
            file.write('p[3] pressure in [Pa] is: %.3g\n' %pxn[3])
            file.write('p[-1] pressure in [Pa] is: %.3g\n' %pxn[-1])
            file.write('p[-2] pressure in [Pa] is: %.3g\n' %pxn[-2])
            file.write('p[-3] pressure in [Pa] is: %.3g\n' %pxn[-3])
            file.write('max pressure is:  %.3g\n' %numpy.max(pxn))
            file.write('min pressure is:  %.3g\n' %numpy.min(pxn))
            file.write('momentum cfl number is:  %.3g\n' %cflmn)
            file.write('energy cfl number is:  %.3g\n' %cflen)
            file.write('anode temperature addition [K]:  %.4g\n' %Tpa)
            file.write('cathode temperature addition [K]:  %.4g\n' %Tpc)
            file.write('---------------------------------------\n')
            file.close()
	#print graphs
        if iterations in reporting_interval2:
            plot( xc,numpy.round(Tn[1:-1],decimals=10),\
                 'grid location [m]','Temperature [K]',  'time= %.4g sec' %time)
            pyplot.savefig('./plots/%.5g_temperature.png' %iterations)
            #
            plot( xb,numpy.round(uxn[1:-1],decimals=5),\
                 'grid location [m]','velocity [m/s]',  'time= %.4g sec' %time)
            pyplot.savefig('./plots/%.5g_velocity.png' %iterations)
            #
            plot( xc,numpy.round(rhoxn[1:-1],decimals=4),\
                 'grid location [m]','density [kg/m3]',  'time= %.4g sec' %time)
            pyplot.savefig('./plots/%.5g_density.png'%iterations)
            #
            plot( xc,numpy.round(pxn[1:-1],decimals=0),\
                 'grid location [m]','pressure [Pa]',  'time= %.4g sec' %time)
            pyplot.savefig('./plots/%.5g_pressure.png'%iterations)
            #
            plot( xc,numpy.round(nexn[1:-1],decimals=0),\
                 'grid location [m]','ne [1/m3]',  'time= %.4g sec' %time)
            pyplot.savefig('./plots/%.5g_ne.png'%iterations)
            #
            plot( xc,numpy.round(econdxn[1:-1],decimals=0),\
                 'grid location [m]','econd [S/m]',  'time= %.4g sec' %time)
            pyplot.savefig('./plots/%.5g_econdavg.png'%iterations)
            #
            plot( xc,numpy.round(rhoB[1:-1],decimals=4),\
                 'grid location [m]','Boron Density [kg/m3]','time= %.4g sec' %time)
            pyplot.savefig('./plots/%.5g_rhoB.png' %iterations)
            #
            plot( xc,numpy.round(rhoNi[1:-1],decimals=4),\
                 'grid location [m]','Nickel Density [kg/m3]','time= %.4g sec' %time)
            pyplot.savefig('./plots/%.5g_rhoNi.png' %iterations)
            #
            plot( xc,numpy.round(rhoCo[1:-1],decimals=4),\
                 'grid location [m]','Cobalt Density [kg/m3]','time= %.4g sec' %time)
            pyplot.savefig('./plots/%.5g_rhoCo.png' %iterations)
            #
            plot( xc,numpy.round(rhoN[1:-1],decimals=4),\
                 'grid location [m]','Nitrogen Density [kg/m3]','time= %.4g sec' %time)
            pyplot.savefig('./plots/%.5g_rhoN.png' %iterations)
	#check for anode steady state:
       # if numpy.abs(Tnref[0]-Tnref_old[0]) <l2_target and iterations>1e6:
           # Tnref[0] = Tnref_old[0]
           # print('anode has reached steady state')
	#check for overall steady state
        if (numpy.all(numpy.abs(Tnref[:]-Tnref_old[:]) < l2_target)==True):
            plot( xc,numpy.round(Tn[1:-1],decimals=10),\
                 'grid location [m]','Temperature [K]',  'time= %.4g sec' %time)
            pyplot.savefig('./plots/%.5g_temperature.png' %iterations)
            #
            plot( xb,numpy.round(uxn[1:-1],decimals=5),\
                 'grid location [m]','velocity [m/s]',  'time= %.4g sec' %time)
            pyplot.savefig('./plots/%.5g_velocity.png' %iterations)
            #
            plot( xc,numpy.round(rhoxn[1:-1],decimals=4),\
                 'grid location [m]','density [kg/m3]',  'time= %.4g sec' %time)
            pyplot.savefig('./plots/%.5g_density.png'%iterations)
            #
            plot( xc,numpy.round(pxn[1:-1],decimals=0),\
                 'grid location [m]','pressure [Pa]',  'time= %.4g sec' %time)
            pyplot.savefig('./plots/%.5g_pressure.png'%iterations)
            #
            plot( xc,numpy.round(nexn[1:-1],decimals=0),\
                 'grid location [m]','ne [1/m3]',  'time= %.4g sec' %time)
            pyplot.savefig('./plots/%.5g_ne.png'%iterations)
            #
            plot( xc,numpy.round(econdxn[1:-1],decimals=0),\
                 'grid location [m]','econd [S/m]',  'time= %.4g sec' %time)
            pyplot.savefig('./plots/%.5g_econdavg.png'%iterations)
            file.write('entire domain has reached steady state\n')
            file.write('terminal time is %.4g\n' %time)
            file.write('end of computation, see you later!\n')
            break
        #check for end of computation
        if (time >= t_terminal):
            file.write('max computational time reached.if steady state not achieved, increase terminal time\n')
    return iterations,time	
