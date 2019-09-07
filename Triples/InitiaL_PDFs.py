import numpy as np 
import scipy.stats as st
import math 
from scipy.integrate import simps
import datetime
import sys

#Making sure the numbers in the text file outputterd are not truncated
np.set_printoptions(precision=5)
np.set_printoptions(threshold=sys.maxsize)

#datetime module is used to stamp the simulation files
now = datetime.datetime.now()
now = now.strftime("%d-%m-%Y_%H-%M-%S")

'''Defining constants'''
#All masses are defined in solar masses 
#Take M1 To be 2 solar masses 
min_star_mass = 0.08
M1= 2 
Rsun = 695510*10**3
Msun = 2*10**30
G = 6.67*10**(-11)


#This is the general method for defining own distribution
'''Block of code used to generate thermal PDF for eccentricities'''
class thermal(st.rv_continuous):
    def _pdf(self,x):
        return p(x)
def p(x):
    #Assume f(e)= 2e distribution for eccentricities
    return 2*x
def normalisation(x):
    return simps(p(x),x)
    


def uniform_dist(start,end,num_points):
    '''Generates uniform PDF'''
    random_draw = np.random.uniform(start,end,num_points)
    return(random_draw)
    
def normal_dist(num_points,mu,sigma):
    '''Generates Normal PDF'''
    random_draw = np.random.normal(mu,sigma,num_points)
    return(random_draw)

def T14_orbital_draw(num_samples,mean,std_dev):
    '''Draws orbital seperations from T14 distribution'''
    #NOTE: T14 orbital distribution is defined in log(Period) as mu = 5 std = 2.3
    P_samples = normal_dist(num_samples,mean,std_dev)
    #Convert period distribution to orbital seperation dsitribution
    P_samples = (10**P_samples)*24*60*60
    a_samples = (G*M1*Msun/(4*math.pi**2.) * P_samples**2.)**(1./3.)
    a_normal = np.log10(a_samples/Rsun)
    return(a_normal)

def OBin_orbital_draw(num_samples,lower_bound,upper_bound):
    '''Draws orbital seperations from OBin distribution'''
    a_uniform = uniform_dist(lower_bound,upper_bound,num_samples)
    return a_uniform

def Sort_Orbital_Seperations(samples):
    '''Sorts array of orbital seperations into triple configuation of ain/aout'''
    ain = []
    aout  = []
    while len(samples)!=0:
        draw1 = samples[0]
        draw2 = samples[1]
        samples = np.delete(samples,0)
        samples = np.delete(samples,0)
        if draw1>draw2 and draw1 != draw2:
            aout.append(draw1)
            ain.append(draw2)
        elif draw2>draw1 and draw1 != draw2:
            aout.append(draw2)
            ain.append(draw1)
    
    return(np.array(ain),np.array(aout))

def uniform_mass_draw(num_points):
    '''Draws mass ratios from uniform mass disttribution'''
    qoutlist = uniform_dist(min_star_mass/(2*M1) ,1,num_points)
    return qoutlist

def oneoverq_mass_draw(num_points):
    '''Draws mass ratios from 1/q mass distribution'''
    #We generate uniform distribtuion in log(q), this generates a 1/q distribution in q
    qoutlist = uniform_dist(math.log10(min_star_mass/(2*M1)),math.log10(1),num_points)
    qoutlist = 10**qoutlist
    return qoutlist

def eccentricity_thermal_draw(num_points):
    '''Defines thermal distribution and draws eccentricities from it'''
    #Defining PDF
    eccentricity_dist= thermal(name="Eccentricity",a=0.0)
    eccentricities = np.linspace(0,1,1000)
    scale_factor = 1/normalisation(eccentricities)
    samples_eccentricities = eccentricity_dist.rvs(size = num_points,scale = scale_factor)
    return samples_eccentricities

def inclinations_draw(num_points):
    '''Draws inclinations from distribution uniform in cos(i) (between pi and -pi)'''
    inclination_dist = uniform_dist(1,-1,num_points)
    inclinations = np.arccos(inclination_dist)
    return(inclinations)
    
def Save_Data(filename,ain,aout,qoutlist,eccentricites,inclinations):
    with open(filename,'a',newline='') as output:
        output.write('Format: ain, aout,mass ratios,eccentricities,inclinations \n')
        output.write('#'+str(10.00**ain)+'\n')
        output.write('#'+str(10.00**aout)+'\n')
        output.write('#'+str(qoutlist)+'\n') 
        output.write('#'+str(eccentricites)+'\n')
        output.write('#'+str(inclinations)+'\n')
        


#samples = T14_orbital_draw(50000,5,2.3)
#ain,aout = Sort_Orbital_Seperations(samples)
#qoutlist = uniform_mass_draw(100)
#eccentricity_draws = eccentricity_thermal_draw(100)
#inclination_draws = inclinations_draw(100)
#Save_Data('Triple_Simulation_Data'+str(now)+'.txt',ain,aout,qoutlist,eccentricity_draws,inclination_draws)
#














    

                    



    
    
    