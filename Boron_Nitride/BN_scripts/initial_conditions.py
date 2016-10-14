#This script initializes arrays for use
#Luis Martinez 
#September 29,2016
#import libraries
import numpy
from parameters import R_N,k_an,cp_N,rho_N,Lgap,nx_gap,I_arc,Ran
from electrical_properties import get_electric_potential
from check_interval import get_report_interval
#
'''
This script contains the following functions

initial_conditions.get_domain
initial_conditions.get_initial_conditions
initial_conditions.get_initial_concentrations
electrical_properties.get_electric_potential

'''
#get domain:
def get_domain(Ld, n_gap):
    ''' This function obtains the domain based on a 
    staggered grid. The domain is divided into three parts:
    (1) anode region ghost cell, (2) gap, (3) cathode region ghost cell.
    
    Parameters:
    ----------
    Ld: Domain Length
    n_gap: number of boundary nodes in gap
    
    Returns:
    --------
    delx_gap: x-direction spacing
    xb_gap: cell boundary nodes
    xc_gap: cell center nodes   
    '''
    #staggered grid approach
    #gap region
    xl_gap = 0.0
    xr_gap = Ld
    delx_gap = (xr_gap-xl_gap)/(n_gap-1.0)
    xb_gap = numpy.linspace(xl_gap, xr_gap, n_gap)
    xc_gap = xb_gap[:-1]+(delx_gap*0.5)
    #
    return delx_gap, xb_gap, xc_gap
#----------------------------------------------------
dx,xb,xc = get_domain(Lgap,nx_gap)
#----------------------------------------------------
#
def get_initial_conditions(x1,x2,I_input,R_anode):
	'''
	- Defines initial arrays
    	- Adds two ghost cells to represent nodes on anode and cathode,
        for boundary conditions
    
    	Parameters:
    	-----------
    	x1: boundary x locations
    	x2: center cell x locations
    	I_input: arc current input
    	R_anode: anode radius
    
    	Returns:
    	-------
    	Arrays!
    	Tx: temperature 
    	jx: current density
    	phix: electric potential
    	cpx: specific heat
    	visc: viscosity (mu)
    	ux: velocity
    	kx: thermal conductivity
    	rhox: density
    	Px: pressure
    	hx: enthalpy
    	econdx: electrical conductivity
    	nex: number density of electrons
    
   	 Consistent with the staggered grid approach, the velocity is calculated on the
    	cell boundaries and all of the other values are computed on the cell centers.
	'''
	Tref = 300.
	a = numpy.size(x1)+2 #adds 2 ghost cells to array
	b = numpy.size(x2)+2 #adds 2 ghost cells to array
	#intialize arrays
	ux = numpy.zeros((a),dtype=float) #velocity
	jx = numpy.zeros((b),dtype=float) #current density
	phix = numpy.zeros((b),dtype=float) #electric potential
	econdx = numpy.zeros((b),dtype=float) #electrical conductivity
	rhox = numpy.zeros((b),dtype=float) #density
	nex = numpy.zeros((b),dtype=float) #electron density
	cpx = numpy.zeros((b),dtype=float) #specific heat
	Tx = numpy.zeros((b),dtype=float) #temperature
	mux = numpy.zeros((b),dtype=float) #shear viscosity
	kx = numpy.zeros((b),dtype=float) #thermal conductivity
	px = numpy.zeros((b),dtype=float) #pressure
	hx = numpy.zeros((b),dtype=float) #enthalpy
	Rsx = numpy.zeros((b),dtype=float) #gas constant
	#initialize values
	jx[:] = I_input/(numpy.pi*(R_anode**2))
	econdx[:] = I_input
	rhox[:] = rho_N
	nex[:] = 1e6
	cpx[:] = cp_N
	Tx[:] = Tref
	mux[:] = 1.7e-6 #nitrogen at room temperature
	kx[:] = 0.02583
	kx[0] = k_an
	px[:] = rhox[:]*Tx[:]*R_N
	hx[:] = cpx[:]*Tx[:]
	#return
	return Tx, jx, phix, cpx, mux, ux, kx, rhox, px, hx, econdx, nex
#-------------------------------------------------------------------------
Tn,jn,phin,cpn,mun,uxn,kn,\
rhoxn,pxn,hn,econdxn,nexn = get_initial_conditions(xb,xc,I_arc,Ran)
#------------------------------------------------------------------------
numx = numpy.size(pxn)
#------------------------------------------------------------------------
def initial_concentrations(density):
    '''
    Defines initial concentrations of all species in the domain

    Parameters:
    -----------
    density

    Returns:
    --------
    cs1: Boron
    cs2: Nickel
    cs3: Cobalt
    cg1: Nitrogen
    '''
    cs1 = numpy.zeros_like(density)
    cs2 = numpy.zeros_like(density)
    cs3 = numpy.zeros_like(density)
    cg1 = numpy.zeros_like(density)
    #
    cs1[:] = 1e-14
    cs2[:] = 1e-14
    cs3[:] = 1e-14
    cg1[:] = 1. - cs1[:] - cs2[:] - cs3[:]
    #
    return cs1,cs2,cs3,cg1
#-----------------------------------------------------------------
cB,cNi,cCo,cN = initial_concentrations(rhoxn)
#-----------------------------------------------------------------
phin = get_electric_potential(jn.copy(),econdxn.copy(), phin.copy(), xc.copy(), dx)
#-----------------------------------------------------------------
max_iteration = 1e10
step_interval = 5e3
step_interval2 = 10e3
check_iterations_size = int(max_iteration/step_interval)
report_int1 = numpy.zeros((check_iterations_size))
report_int2 = numpy.zeros((check_iterations_size))
rint1,rint2  = get_report_interval(report_int1,report_int2,max_iteration,step_interval,step_interval2)
