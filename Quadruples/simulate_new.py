import numpy as np 
import datetime
import sys
import logging
import multiprocessing
import simulationQuad


logging.basicConfig(level=logging.DEBUG,
                    format='(%(threadName)-10s) %(message)s',
                    )
np.set_printoptions(precision=5)
np.set_printoptions(threshold=sys.maxsize)
#Datetime module is used for timestamps
now = datetime.datetime.now()
now = now.strftime("%d-%m-%Y_%H-%M-%S")



def load_data(filename,mass_ratio_filename):
    '''Loading in orbital parameters drawn from defined PDFs'''
    with open(filename,'r') as file:
        contents = file.read()
        
    ain = contents.replace('\n','').split('#')[1]
    amid = contents.replace('\n','').split('#')[2]
    aout = contents.replace('\n','').split('#')[3]
    inclinations = contents.replace('\n','').split('#')[4]
    samples_eccentricities = contents.replace('\n','').split('#')[5]
    q_shape = contents.replace('\n','').split('#')[6]

    ain = np.array([float(i) for i in ain.replace('[','').replace(']','').split()])
    amid = np.array([float(i) for i in amid.replace('[','').replace(']','').split()])
    aout = np.array([float(i) for i in aout.replace('[','').replace(']','').split()])
    inclinations = np.array([float(i) for i in inclinations.replace('[','').replace(']','').split()])
    samples_eccentricities = np.array([float(i) for i in samples_eccentricities.replace('[','').replace(']','').split()])
    q_shape  = [float(i) for i in q_shape.replace('(','').replace(')','').split(',')]
                                
    q_final = np.loadtxt(mass_ratio_filename).reshape((int(q_shape[0]),int(q_shape[1])))

    return(ain,amid,aout,inclinations,samples_eccentricities,q_final,q_shape)

def calculate(i):
    '''Computes stability for one section of the parameter space (split by a_ratios)'''
    global q_final_counter
    result = simulationQuad.simulate(i=i+1,a_ratio1=a_ratio1[i],a_ratio2=a_ratio2[i],qoutlist=q_final,samples_eccentricities = samples_eccentricities,
			      inclinations = inclinations)
    a_ratio_counter[i] = result[0]
    for i in range(len(result[1])):
        primary_eccentricity_counter[i] += result[1][i]
    for i in range(len(result[2])):
        secondary_eccentricity_counter[i] += result[2][i]
    q_final_counter += result[3]
    return(list(result[1]))

def save_data(filename,massresult_filename):
    '''Output counter variables from simulation into text file'''

    with open(filename,'a') as output:
        output.write(str(np.array(a_ratio_counter[:],dtype='int64'))+'\n')
        output.write('#'+str(np.array(primary_eccentricity_counter[:],dtype = 'int64'))+'\n')
        output.write('#'+str(np.array(secondary_eccentricity_counter[:],dtype = 'int64'))+'\n')
    
    np.savetxt(massresult_filename,q_final_counter)



ain,amid,aout,inclinations,samples_eccentricities,q_final,q_shape = load_data('Initial_data.txt','Initial_mass_ratios.txt')



a_ratio1 = amid/ain
a_ratio2  = aout/amid

'''Set up final counter arrays as multiprocessing arrays which all processes can access'''
a_ratio_counter = multiprocessing.Array("d", [0]*len(ain))
primary_eccentricity_counter = multiprocessing.Array("d", [0]*len(samples_eccentricities))
secondary_eccentricity_counter = multiprocessing.Array("d", [0]*len(samples_eccentricities))
q_final_counter = multiprocessing.Array("d",int(q_shape[0]*q_shape[1]))
q_final_counter = np.frombuffer(q_final_counter.get_obj())
q_final_counter = q_final_counter.reshape((int(q_shape[0]),int(q_shape[1])))


'''Deals with the parallelisation'''
#If want to use smaller number of cores, change multiprocessing.cpu_count() to "processes = n" 
#where n is the number of cores to be used
#By default set to the maximum number of cores
pool = multiprocessing.Pool(multiprocessing.cpu_count())
#Assign each process to use parameter space occupied by one element of a_ratio 
a_ratio_index = [i for i in range(len(a_ratio1))] 
closest_counter = pool.map(calculate,a_ratio_index)

q_final_counter = np.array(q_final_counter,dtype='int64')
save_data('Simulation_results.txt','mass_array_results_counter.txt')




