# This script contains all of the constant variables in the domain
#Luis Martinez
#September 29, 2016 
'''
This script contains all of the constants used!
'''
import scipy.constants
import numpy
# Input current (Amperes)
I_arc = 60.

#Collision Radii, set equal to Van Der Waal radii
RSCD_N_B = 1.55+1.92 # Boron-Nitrogen
RSCD_N_Ni = 1.55 + 1.63 # Nitrogen-Nickel (these last two are approximated values)
RSCD_N_Co = 1.55 + 1.63 # Nitrogen-Cobalt 

#physical constants
Ckb = scipy.constants.k  #Boltzmann Constant J/K
Cme = scipy.constants.m_e #Elementary Electron Mass kg
Ch = scipy.constants.h; # Planck's constant Js
CNA = scipy.constants.N_A #avogadros number
Ce = scipy.constants.e #elementary charge, C
CA0 = 1.20173e6 #A/[m2 K2] #Richardson constant
CAG = 0.5*CA0 #factor for thermionic emission at cathode end

#Boron
rho_B = 2460. #kg/m3
M_B = 10.811e-3/CNA #kg, mass of 1 atom of pure boron (also = Ckb/M)
mB = 10.811 # molar mass g/mol
Ips1 = 800.6e3/CNA #J, energy of first ionization
cp_B = 1030. # J/(Kg K)
R_B = 8314.0/mB #specific gas constant
Boron_workfunction = 4.45 #eV
dH_B = 507e6/mB #heat of vaporization J/Kg
Uiz_B = 8.2980 #eV
k_B = 27. #W/(mK)

#Nickel
rho_Ni = 8908
M_Ni = 58.6934e-3/CNA #kg, mass of 1 atom of pure nickel (also = Ckb/M)
mNi = 58.6934 # molar mass g/mol
Ips2 = 737.1e3/CNA #J, energy of first ionization
cp_Ni = 445. # J/(Kg K)
R_Ni = 8314.0/mNi #specific gas constant
Nickel_workfunction = 5.01 #eV
dH_Ni = 378e6/mNi #heat of vaporization J/Kg
Uiz_Ni = 7.6398 #eV
k_Ni = 91. #W/(m K)

#Cobalt
rho_Co = 8900.
M_Co = 58.9332/CNA #kg, mass of 1 atom of pure cobalt (also = Ckb/M)
mCo = 58.9332 # molar mass g/mol
Ips3 = 760.4e3/CNA #J, energy of first ionization
cp_Co = 418.68# J/(Kg K)
R_Co = 8314.0/mCo #specific gas constant
Cobalt_workfunction = 5.0 #eV
dH_Co = 375e6/mCo #heat of vaporization J/Kg
Uiz_Co = 7.8810 #eV
k_Co = 100. #W/(m K)

#Nitrogen Gas
rho_N = 1.25
M_N = (2*14.0067)/CNA # kg, mass of 1 atom of Helium
mN = (2*14.0067) #molar mass, g/mol
Ipg1 = 1402.3e3/CNA; # J, energy of first ionization
R_N = 8314.0/mN #specific gas constant
cp_N = 1039. #J/Kg-K
Uiz_N = 14.5341 #eV
k_N = 0.02583 #W/(m K)

#assume composition of ingot (grams)
WB = 2.00
WNi = 0.5*(WB*(1-0.99))/0.99
WCo = 0.5*(WB*(1-0.99))/0.99
Wtotal = WB+WNi+WCo

#mass fractions
wB = WB/(Wtotal)
wNi = WNi/(Wtotal)
wCo = WCo/(Wtotal)

#molar fractions
xtotal = WB/mB + WNi/mNi + WCo/mCo
xB = (WB/mB)/xtotal
xNi = (WNi/mNi)/xtotal
xCo = (WCo/mCo)/xtotal
xN = 0.

#molar mass of evaporated anode material
molmass_evap = 1/((wB/mB)+(wNi/mNi)+(wCo/mCo))

#anode 
rho_an_max = 0.99*(rho_B) + 0.005*(rho_Ni+rho_Co)
Rs_an = 0.99*(R_B) + 0.005*(R_Ni+R_Co)
MW_evap = 0.99*M_B  + 0.005*(M_Ni + M_Co)
k_an = 0.99*k_B  + 0.005*(k_Ni + k_Co)
cp_an = 0.99*cp_B  + 0.005*(cp_Ni + cp_Co)
Uiz_an = 0.99*Uiz_B  + 0.005*(Uiz_Ni + Uiz_Co)
dH_an = 0.99*dH_B  + 0.005*(dH_Ni + dH_Co)
anode_workfunction = 0.99*Boron_workfunction +\
                    0.005*(Nickel_workfunction+Cobalt_workfunction)

#cathode: Tungsten
cathode_workfunction = 4.5 #eV
k_cath = 173. # W/(m K) 
rho_cath_max = 19250. #kg/m3 
cp_cath = 1340. #J/(kg K)
Uiz_cath = 7.8640 #eV
 
#Domain information
Lgap = 0.004 #interelectrode gap [meters]
nx_gap = 121 #number of total boundaries
#Rc = (12.5/2.) * 1./1000. # cathode radius [meters]
Ran = 6.35/2. * 1./1000.  #anode radius [meters] 
l2_target = 1e-7 #error tolerance for pressure equation
dt = 4e-9 #{sec} initial time step
#
#array for l2 norm of error in pressure equation
#done here to avoid using numpy with numba
l2_err = numpy.zeros(20000)
