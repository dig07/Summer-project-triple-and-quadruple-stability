import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import math
import matplotlib
resolution=10

plt.rc('legend',**{'fontsize':22})
matplotlib.rcParams['xtick.labelsize'] = 22
matplotlib.rcParams['ytick.labelsize'] = 22 
matplotlib.rcParams['axes.labelsize'] = 28


#Function used to plot the histograms 
def plot_piecewise_func(boundaries,values):
    x = []
    y = []
    values = list(values)
    for section in boundaries:
        xrange = list(np.linspace(section[0],section[1],resolution))
        yrange = [values[boundaries.index(section)]]*len(xrange)
        x.extend(xrange)
        y.extend(yrange)
    return(x,y)
    
#Load the simulation details file
with open('simulation_details_QUADRUPLE_uniform05-08-2019_22-20-51.txt','r') as file:
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

#Load in mass ratio array
q_final = np.loadtxt('Q_RATIO_ARRAY_details_QUADRUPLE_uniform05-08-2019_22-20-51.txt').reshape((int(q_shape[0]),int(q_shape[1])))


#Load the simulation results 
resulttext = 'Test_STABLE_Qratio_uniform_uniform05-08-2019_22-23-16.txt'
with open(resulttext,'r')as file:
    results = file.read()

a_ratio= results.replace('\n','').split('#')[0]
primary_eccentricity_counter = results.replace('\n','').split('#')[1]
secondary_eccentricity_counter = results.replace('\n','').split('#')[2]
                                                
a_ratio=[float(i) for i in a_ratio.replace('[','').replace(']','').split()]
primary_eccentricity_counter=[float(i) for i in primary_eccentricity_counter.replace('[','').replace(']','').split()]
secondary_eccentricity_counter = np.array([float(i) for i in secondary_eccentricity_counter.replace('[','').replace(']','').split()])

#load in mass ratio counter array
q_final_counter = np.loadtxt('Q_RATIO_ARRAY_STABLE_QUADRUPLE_uniform05-08-2019_22-23-16.txt').reshape((int(q_shape[0]),int(q_shape[1])))

#Sort primary and secondary mass ratios and counts into seperate arrays
primary_mass_ratios = []
primary_mass_ratio_counts = []
secondary_mass_ratios = []
secondary_mass_ratio_counts = []

for i in range(len(q_final_counter[0,:])):
    primary_mass_ratios.append(q_final[0,:][i])
    primary_mass_ratio_counts.append(q_final_counter[0,:][i])
    secondary_mass_ratios.extend(list(q_final[1:,i]))
    secondary_mass_ratio_counts.extend(list(q_final_counter[1:,i]))
    

primary_mass_ratios = np.array(primary_mass_ratios)
primary_mass_ratio_counts = np.array(primary_mass_ratio_counts)
secondary_mass_ratios = np.array(secondary_mass_ratios)
secondary_mass_ratio_counts = np.array(secondary_mass_ratio_counts)



     

    


'''For plotting secondary mass distributions'''

df = pd.DataFrame({'mass ratios considered':primary_mass_ratios,'mass ratio counts':primary_mass_ratio_counts})
df2 = pd.DataFrame({'mass ratios considered':secondary_mass_ratios,'mass ratio counts':secondary_mass_ratio_counts})

count, bins, ignored = plt.hist(primary_mass_ratios, 30, normed=True,visible=False)
bins = list(bins)

#Sections for piecewise function 
boundaries = []
for boundary in bins:
    if boundary == bins[len(bins)-1]:
        pass
    else:
        boundaries.append((boundary,bins[(bins.index(boundary)+1)]))
        
    

#Scale initial counts down by how many stable orbits 
secondary_counts = []
for section in boundaries:
    subsection = df[(df['mass ratios considered'] > section[0]) & (df['mass ratios considered'] < section[1])]
    num_unique_orbits = len(subsection['mass ratios considered'])*len(a_ratio)*q_shape[0]*(len(inclinations)**2)*(len(samples_eccentricities)**2)
    num_stable_orbits = np.sum(subsection['mass ratio counts'])
    if num_unique_orbits == 0 and num_stable_orbits == 0:
        secondary_counts.append(0)
    else:
        factor = num_stable_orbits/num_unique_orbits
        index = boundaries.index(section)
        secondary_counts.append(count[index]*factor)

plt.figure(figsize=(20,10))

#Accounting for the vertical lines 
xoriginal,yoriginal=plot_piecewise_func(boundaries,count)
plt.plot(xoriginal,yoriginal,'b--',label="mass ratios (parent)")
xstable,ystable = plot_piecewise_func(boundaries,secondary_counts)

#Plotting edges (Bit of a bodge)
xmin = boundaries[0][0]
ymin = secondary_counts[0]
xmax = boundaries[(len(boundaries)-1)][1]
ymax = secondary_counts[(len(count)-1)]

plt.vlines(xmin,0,ymin,color='b',linestyle='solid')
plt.vlines(xmax,0,ymax,color='b',linestyle='solid')



ymin = count[0]
ymax = count[(len(count)-1)]

plt.vlines(xmin,0,ymin,color='b',linestyle='--')
plt.vlines(xmax,0,ymax,color='b',linestyle='--')




plt.plot(xstable,ystable,'b',label = "mass ratios (stable)")

plt.legend()
plt.xlabel(r'$\frac{M3}{M1+M2}$')
plt.savefig('Primary mass ratios')

#Secondary mass ratios
count, bins, ignored = plt.hist(secondary_mass_ratios, 30, normed=True,visible=False)
bins = list(bins)

#Sections for piecewise function 
boundaries = []
for boundary in bins:
    if boundary == bins[len(bins)-1]:
        pass
    else:
        boundaries.append((boundary,bins[(bins.index(boundary)+1)]))
        
    

#Scale initial counts down by how many stable orbits 
secondary_counts = []
for section in boundaries:
    subsection = df2[(df2['mass ratios considered'] > section[0]) & (df2['mass ratios considered'] < section[1])]
    num_unique_orbits = len(subsection['mass ratios considered'])*len(a_ratio)*(len(inclinations)**2)*(len(samples_eccentricities)**2)
    num_stable_orbits = np.sum(subsection['mass ratio counts'])
    if num_unique_orbits == 0 and num_stable_orbits == 0:
        secondary_counts.append(0)
    else:
        factor = num_stable_orbits/num_unique_orbits
        index = boundaries.index(section)
        secondary_counts.append(count[index]*factor)

plt.figure(figsize=(20,10))

#Accounting for the vertical lines 
xoriginal,yoriginal=plot_piecewise_func(boundaries,count)
plt.plot(xoriginal,yoriginal,'b--',label="mass ratios (parent)")
xstable,ystable = plot_piecewise_func(boundaries,secondary_counts)

#Plotting edges (Bit of a bodge)
xmin = boundaries[0][0]
ymin = secondary_counts[0]
xmax = boundaries[(len(boundaries)-1)][1]
ymax = secondary_counts[(len(count)-1)]

plt.vlines(xmin,0,ymin,color='b',linestyle='solid')
plt.vlines(xmax,0,ymax,color='b',linestyle='solid')



ymin = count[0]
ymax = count[(len(count)-1)]

plt.vlines(xmin,0,ymin,color='b',linestyle='--')
plt.vlines(xmax,0,ymax,color='b',linestyle='--')




plt.plot(xstable,ystable,'b',label = "mass ratios (stable)")

plt.legend()
plt.xlabel(r'$\frac{M4}{M1+M2+M3}$')
plt.savefig('Secondary mass ratios')


     
        

        



















'''For plotting the eccentricity parent and stable distributions'''

df = pd.DataFrame({'eccentricities considered':samples_eccentricities,'eccentricity counts':primary_eccentricity_counter})

count, bins, ignored = plt.hist(samples_eccentricities, 30, normed=True,visible=False)
bins = list(bins)

#Sections for piecewise function 
boundaries = []
for boundary in bins:
    if boundary == bins[len(bins)-1]:
        pass
    else:
        boundaries.append((boundary,bins[(bins.index(boundary)+1)]))
        
    

#Scale initial counts down by how many stable orbits 
secondary_counts = []
for section in boundaries:
    subsection = df[(df['eccentricities considered'] > section[0]) & (df['eccentricities considered'] < section[1])]
    num_unique_orbits = len(subsection['eccentricities considered'])*len(a_ratio)*q_shape[0]*q_shape[1]*(len(inclinations)**2)*len(samples_eccentricities)
    num_stable_orbits = np.sum(subsection['eccentricity counts'])
    if num_unique_orbits == 0 and num_stable_orbits == 0:
        secondary_counts.append(0)
    else:
        factor = num_stable_orbits/num_unique_orbits
        index = boundaries.index(section)
        secondary_counts.append(count[index]*factor)

plt.figure(figsize=(20,10))

#Accounting for the vertical lines 
xoriginal,yoriginal=plot_piecewise_func(boundaries,count)
plt.plot(xoriginal,yoriginal,'b--',label="eccentricity (parent)")
xstable,ystable = plot_piecewise_func(boundaries,secondary_counts)

#Plotting edges (Bit of a bodge)
xmin = boundaries[0][0]
ymin = secondary_counts[0]
xmax = boundaries[(len(boundaries)-1)][1]
ymax = secondary_counts[(len(count)-1)]

plt.vlines(xmin,0,ymin,color='b',linestyle='solid')
plt.vlines(xmax,0,ymax,color='b',linestyle='solid')



ymin = count[0]
ymax = count[(len(count)-1)]

plt.vlines(xmin,0,ymin,color='b',linestyle='--')
plt.vlines(xmax,0,ymax,color='b',linestyle='--')




plt.plot(xstable,ystable,'b',label = "eccentricity (stable)")

plt.legend()
plt.xlabel('e')
plt.savefig('Primary Eccentricity dist')






'''For plotting the secondary eccentricity parent and stable distributions'''

df = pd.DataFrame({'eccentricities considered':samples_eccentricities,'eccentricity counts':secondary_eccentricity_counter})

count, bins, ignored = plt.hist(samples_eccentricities, 30, normed=True,visible=False)
bins = list(bins)

#Sections for piecewise function 
boundaries = []
for boundary in bins:
    if boundary == bins[len(bins)-1]:
        pass
    else:
        boundaries.append((boundary,bins[(bins.index(boundary)+1)]))
        
    

#Scale initial counts down by how many stable orbits 
secondary_counts = []
for section in boundaries:
    subsection = df[(df['eccentricities considered'] > section[0]) & (df['eccentricities considered'] < section[1])]
    num_unique_orbits = len(subsection['eccentricities considered'])*len(a_ratio)*q_shape[0]*q_shape[1]*(len(inclinations)**2)*len(samples_eccentricities)
    num_stable_orbits = np.sum(subsection['eccentricity counts'])
    if num_unique_orbits == 0 and num_stable_orbits == 0:
        secondary_counts.append(0)
    else:
        factor = num_stable_orbits/num_unique_orbits
        index = boundaries.index(section)
        secondary_counts.append(count[index]*factor)

plt.figure(figsize=(20,10))

#Accounting for the vertical lines 
xoriginal,yoriginal=plot_piecewise_func(boundaries,count)
plt.plot(xoriginal,yoriginal,'b--',label="eccentricity (parent)")
xstable,ystable = plot_piecewise_func(boundaries,secondary_counts)

#Plotting edges (Bit of a bodge)
xmin = boundaries[0][0]
ymin = secondary_counts[0]
xmax = boundaries[(len(boundaries)-1)][1]
ymax = secondary_counts[(len(count)-1)]

plt.vlines(xmin,0,ymin,color='b',linestyle='solid')
plt.vlines(xmax,0,ymax,color='b',linestyle='solid')



ymin = count[0]
ymax = count[(len(count)-1)]

plt.vlines(xmin,0,ymin,color='b',linestyle='--')
plt.vlines(xmax,0,ymax,color='b',linestyle='--')




plt.plot(xstable,ystable,'b',label = "eccentricity (stable)")

plt.legend()
plt.xlabel('e')
plt.savefig('Secondary Eccentricity dist')



'''This is for plotting the original histogram with probability density on the y axis 
and log(a/R(sun)) on the x axis'''
#Convert back to logs for plotting

ain_final =[]
for i in ain:
    ain_final.append(math.log10(i))

amid_final = []
for i in amid:
    amid_final.append(math.log10(i))


aout_final = []
for i in aout:
    aout_final.append(math.log10(i))
    
#DONE NORMALLY BECAUSE WE ONLY HAVE MEASUREMENT OF AIN
df = pd.DataFrame({'ain_considered':ain_final,'ain_stable_counts':a_ratio})
df2 = pd.DataFrame({'amid_considered':amid_final,'amid_stable_primary_counts':a_ratio})
df3 = pd.DataFrame({'aout_considered':aout_final,'aout_counter':a_ratio})

#AIN
count, bins, ignored = plt.hist(ain_final, 15, normed=True,visible=False)
bins = list(bins)


#Sections for piecewise function 
boundaries = []
for boundary in bins:
    if boundary == bins[len(bins)-1]:
        pass
    else:
        boundaries.append((boundary,bins[(bins.index(boundary)+1)]))
        
    

#Scale initial counts down by how many stable orbits 
secondary_counts = []
for section in boundaries:
    subsection = df[(df['ain_considered'] > section[0]) & (df['ain_considered'] < section[1])]
    num_unique_orbits = len(subsection['ain_considered'])*q_shape[0]*q_shape[1]*(len(inclinations)**2)*(len(samples_eccentricities)**2)
    num_stable_orbits = np.sum(subsection['ain_stable_counts'])
    if num_unique_orbits == 0 and num_stable_orbits == 0:
        secondary_counts.append(0)
    else:
        factor = num_stable_orbits/num_unique_orbits
        index = boundaries.index(section)
        secondary_counts.append(count[index]*factor)

plt.figure(figsize=(20,10))

#Accounting for the vertical lines 
xoriginal,yoriginal=plot_piecewise_func(boundaries,count)
plt.plot(xoriginal,yoriginal,'b--',label="$a_{in}$ (parent)")
xstable,ystable = plot_piecewise_func(boundaries,secondary_counts)

#Plotting edges (Bit of a bodge)
xmin = boundaries[0][0]
ymin = secondary_counts[0]
xmax = boundaries[(len(boundaries)-1)][1]
ymax = secondary_counts[(len(count)-1)]

plt.vlines(xmin,0,ymin,color='b',linestyle='solid')
plt.vlines(xmax,0,ymax,color='b',linestyle='solid')



ymin = count[0]
ymax = count[(len(count)-1)]

plt.vlines(xmin,0,ymin,color='b',linestyle='--')
plt.vlines(xmax,0,ymax,color='b',linestyle='--')




plt.plot(xstable,ystable,'b',label = "$a_{in}$ (stable)")



###Amid
count, bins, ignored = plt.hist(amid_final, 15, normed=True,visible=False)
bins = list(bins)

boundaries = []
for boundary in bins:
    if boundary == bins[len(bins)-1]:
        pass
    else:
        boundaries.append((boundary,bins[(bins.index(boundary)+1)]))
        
        

secondary_counts = []
for section in boundaries:
    subsection = df2[(df2['amid_considered'] > section[0]) & (df2['amid_considered'] < section[1])]
    num_unique_orbits = len(subsection['amid_considered'])*q_shape[0]*q_shape[1]*(len(inclinations)**2)*(len(samples_eccentricities)**2)
    num_stable_orbits = np.sum(subsection['amid_stable_primary_counts'])
    if num_unique_orbits == 0 and num_stable_orbits == 0:
        secondary_counts.append(0)
    else:
        factor = num_stable_orbits/num_unique_orbits
        index = boundaries.index(section)
        secondary_counts.append(count[index]*factor)

xoriginal,yoriginal=plot_piecewise_func(boundaries,count)
plt.plot(xoriginal,yoriginal,'r--',label="$a_{mid}$ (parent)")
xstable,ystable = plot_piecewise_func(boundaries,secondary_counts)


#Plotting edges (Bit of a bodge)
xmin = boundaries[0][0]
ymin = secondary_counts[0]
xmax = boundaries[(len(boundaries)-1)][1]
ymax = secondary_counts[(len(count)-1)]

plt.vlines(xmin,0,ymin,color='r',linestyle='solid')
plt.vlines(xmax,0,ymax,color='r',linestyle='solid')



ymin = count[0]
ymax = count[(len(count)-1)]

plt.vlines(xmin,0,ymin,color='r',linestyle='--')
plt.vlines(xmax,0,ymax,color='r',linestyle='--')




plt.plot(xstable,ystable,'r',label="$a_{mid}$ (stable)")




###AOUT
count, bins, ignored = plt.hist(aout_final, 15, normed=True,visible=False)
bins = list(bins)

boundaries = []
for boundary in bins:
    if boundary == bins[len(bins)-1]:
        pass
    else:
        boundaries.append((boundary,bins[(bins.index(boundary)+1)]))
        
        

secondary_counts = []
for section in boundaries:
    subsection = df3[(df3['aout_considered'] > section[0]) & (df3['aout_considered'] < section[1])]
    num_unique_orbits = len(subsection['aout_considered'])*q_shape[0]*q_shape[1]*(len(inclinations)**2)*(len(samples_eccentricities)**2)
    num_stable_orbits = np.sum(subsection['aout_counter'])
    if num_unique_orbits == 0 and num_stable_orbits == 0:
        secondary_counts.append(0)
    else:
        factor = num_stable_orbits/num_unique_orbits
        index = boundaries.index(section)
        secondary_counts.append(count[index]*factor)

xoriginal,yoriginal=plot_piecewise_func(boundaries,count)
plt.plot(xoriginal,yoriginal,'g--',label="$a_{out}$ (parent)")
xstable,ystable = plot_piecewise_func(boundaries,secondary_counts)


#Plotting edges (Bit of a bodge)
xmin = boundaries[0][0]
ymin = secondary_counts[0]
xmax = boundaries[(len(boundaries)-1)][1]
ymax = secondary_counts[(len(count)-1)]

plt.vlines(xmin,0,ymin,color='g',linestyle='solid')
plt.vlines(xmax,0,ymax,color='g',linestyle='solid')



ymin = count[0]
ymax = count[(len(count)-1)]

plt.vlines(xmin,0,ymin,color='g',linestyle='--')
plt.vlines(xmax,0,ymax,color='g',linestyle='--')




plt.plot(xstable,ystable,'g',label="$a_{out}$ (stable)")





plt.legend()
plt.xlabel(r'$log(\frac{a}{R_\odot}$)')
plt.savefig('Orbital seperation distribution)


                             

                            
                            
                            
                            
                            
                            
                            
                            
                    