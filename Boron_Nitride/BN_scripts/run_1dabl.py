#this script connects all functions and runs the solution
from parameters import *
from solution import get_solution
from initial_conditions import Tn,jn,phin,cpn,mun,uxn,kn,rhoxn,pxn,hn,econdxn,nexn,cB,cNi,cCo,cN,rint1,rint2,dx
#
itnum,ctime = get_solution(hn,Tn,Tn,pxn,pxn,rhoxn,rhoxn,uxn,cB,cNi,cCo,cN,phin,jn,econdxn,nexn,dx,dt,125,0.75,1e5,1.0,rint1,rint2)
