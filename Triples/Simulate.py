import numpy as np 
import datetime
import math
import sys
import multiprocessing
import simulation

np.set_printoptions(precision=5)
np.set_printoptions(threshold=sys.maxsize)
#Datetime module is used for timestamps
now = datetime.datetime.now()
now = now.strftime("%d-%m-%Y_%H-%M-%S")

def load_data(filename):
    '''Loading in orbital parameters drawn from defined PDFs'''
    with open(filename,'r') as file:
        contents = file.read()
    ain = contents.replace('\n','').split('#')[1]
    aout = contents.replace('\n','').split('#')[2]
    qoutlist = contents.replace('\n','').split('#')[3]
    eccentricities= contents.replace('\n','').split('#')[4]
    inclinations = contents.replace('\n','').split('#')[5]
                                   
    ain = np.array([float(i) for i in ain.replace('[','').replace(']','').split()])
    aout = np.array([float(i) for i in aout.replace('[','').replace(']','').split()])
    qoutlist = np.array([float(i) for i in qoutlist.replace('[','').replace(']','').split()])
    samples_eccentricities = np.array([float(i) for i in eccentricities.replace('[','').replace(']','').split()])
    inclinations = np.array([float(i) for i in inclinations.replace('[','').replace(']','').split()])
    a_ratio = aout/ain
    return(a_ratio,qoutlist,samples_eccentricities,inclinations)

def calculate(i):
    '''Computes stability for one section of the parameter space (split by a_ratios)'''
    result = simulation.simulate(i=i+1,a_ratio=a_ratio[i],qoutlist=qoutlist,samples_eccentricities = samples_eccentricities,
			      inclinations = inclinations)
    #Return and combine results with final counter array
    a_ratio_counter[i] = result[0]
    for i in range(len(result[1])):
        eccentricity_counter[i] += result[1][i]
    for i in range(len(result[2])):
        q_counter[i] += result[2][i]
    #Return eccentricity counter for creation of closest_counter array
    return(list(result[1]))

def save_data(filename):
    '''Output counter variables from simulation into text file'''
    with open(filename,'a') as output:
        output.write(str(np.array(a_ratio_counter[:],dtype='int64'))+'\n')
        output.write('#'+str(np.array(eccentricity_counter[:],dtype = 'int64'))+'\n')
        output.write('#'+str(np.array(q_counter[:],dtype = 'int64'))+'\n')
        output.write('#'+str(np.array(final_closest_counter[:],dtype = 'int64'))+'\n')



a_ratio,qoutlist,samples_eccentricities,inclinations = load_data('Initial_data.txt')

'''Set up final counter arrays as multiprocessing arrays which all processes can access'''
a_ratio_counter = multiprocessing.Array("d", [0]*len(a_ratio))
eccentricity_counter = multiprocessing.Array("d", [0]*len(samples_eccentricities))
q_counter = multiprocessing.Array("d", [0]*len(qoutlist))

'''Deals with the parallelisation'''
#If want to use smaller number of cores, change multiprocessing.cpu_count() to "processes = n" 
#where n is the number of cores to be used
#By default set to the maximum number of cores
pool = multiprocessing.Pool(multiprocessing.cpu_count())
#Assign each process to use parameter space occupied by one element of a_ratio 
a_ratio_index = [i for i in range(len(a_ratio))] 
#Collected eccentricity counter for a(1-e) and a(1-e^2) plot (see report PDF for futher information)
closest_counter = pool.map(calculate,a_ratio_index)

#cleaning up the data for a(1-e) and a(1-e^2) into just one list
final_closest_counter = []
for i in closest_counter:
	final_closest_counter.extend(i)

final_closest_counter = np.array(final_closest_counter, dtype = 'int64')

save_data('Result_STABLE.txt')

