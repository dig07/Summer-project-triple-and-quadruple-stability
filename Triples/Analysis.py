import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import math
import matplotlib
#parameters for plotting
resolution=10
plt.rc('legend',**{'fontsize':22})
matplotlib.rcParams['xtick.labelsize'] = 22
matplotlib.rcParams['ytick.labelsize'] = 22 
matplotlib.rcParams['axes.labelsize'] = 28

#To plot the histogram
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
with open('General_Normal_simulation_details06-08-2019_11-20-11.txt','r') as file:
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



#Used to plot a(1-e) ie the closest point
true_ain = (695510*10**(3))*ain
true_aout = (695510*10**(3))*aout

#Used to plot the closest distance of appriach between the triplt star and inner binary 
total_true_ain_closest = []
total_true_aout_closest  = []

for (ainn,aoutn) in zip(ain,aout):
    li = [ainn*(1+e) for e in samples_eccentricities]
    total_true_ain_closest.extend(li)
    li = [aoutn*(1-e) for e in samples_eccentricities]
    total_true_aout_closest.extend(li)
    
total_true_ain_closest = np.log10(np.array(total_true_ain_closest))
total_true_aout_closest = np.log10(np.array(total_true_aout_closest))


#Load the simulation results 
resulttext = 'Genera_Results_STABLE_normal06-08-2019_11-21-38.txt'
with open(resulttext,'r')as file:
    results = file.read()

a_ratio_counter = results.replace('\n','').split('#')[0]
eccentricity_counter = results.replace('\n','').split('#')[1]
q_final_counter = results.replace('\n','').split('#')[2]
closest_counter= results.replace('\n','').split('#')[3]
                             
#NOTE may need to change the split to ' ' for earlier versions of the code when I was not using arrays
a_ratio_counter=[float(i) for i in a_ratio_counter.replace('[','').replace(']','').split()]
eccentricity_counter = [float(i) for i in eccentricity_counter.replace('[','').replace(']','').split()]
q_final_counter = [float(i) for i in q_final_counter.replace('[','').replace(']','').split()]
closest_counter = [float(i) for i in closest_counter.replace('[','').replace(']','').split()]

#Setting up array for results of a_{in}(1+e_{in})
ain_closest_counter = []
for i in range(len(ain)):
    li = [a_ratio_counter[i]]*len(samples_eccentricities)
    ain_closest_counter.extend(li)
    


'''This is for the plotting the closest approach histogram'''
#total_true_ain_closest = np.array(total_true_ain_closest)
df = pd.DataFrame({'total_closest_ain_considered':total_true_ain_closest,'ain_closest_counts':ain_closest_counter})
df2 = pd.DataFrame({'total_closest_aout_considered':total_true_aout_closest,'aout_closest_counts':closest_counter})

#AIN
count, bins, ignored = plt.hist(total_true_ain_closest, 15, normed=True,visible=False)
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
    subsection = df[(df['total_closest_ain_considered'] > section[0]) & (df['total_closest_ain_considered'] < section[1])]
    num_unique_orbits = len(subsection['total_closest_ain_considered'])*len(qoutlist)*len(samples_eccentricities)*len(inclinations)
    num_stable_orbits = np.sum(subsection['ain_closest_counts'])
    if num_unique_orbits == 0 and num_stable_orbits == 0:
        secondary_counts.append(0)
    else:
        factor = num_stable_orbits/num_unique_orbits
        index = boundaries.index(section)
        secondary_counts.append(count[index]*factor)

plt.figure(figsize=(20,10))

#Accounting for the vertical lines 
xoriginal,yoriginal=plot_piecewise_func(boundaries,count)
plt.plot(xoriginal,yoriginal,'b--',label="$a_{in}(1+e_{in})$ (parent)")
xstable,ystable = plot_piecewise_func(boundaries,secondary_counts)

xfill = [i for i in xstable if i<0.69897]
yfill = [max(count)]*len(xfill)

plt.fill_between(xfill,0,yfill,facecolor='0.5')
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



plt.plot(xstable,ystable,'b',label = "$a_{in}(1+e_{in})$ (stable)")




###AOUT
count, bins, ignored = plt.hist(total_true_aout_closest, 15, normed=True,visible=False)
bins = list(bins)

boundaries = []
for boundary in bins:
    if boundary == bins[len(bins)-1]:
        pass
    else:
        boundaries.append((boundary,bins[(bins.index(boundary)+1)]))
        
        

secondary_counts = []
for section in boundaries:
    subsection = df2[(df2['total_closest_aout_considered'] > section[0]) & (df2['total_closest_aout_considered'] < section[1])]
    #CHECK ASSUMPTION !!!!!!!!!NEED TO CHECK THIS !!!!!--!
    num_unique_orbits = len(subsection['total_closest_aout_considered'])*len(inclinations)*len(qoutlist)
    num_stable_orbits = np.sum(subsection['aout_closest_counts'])
#    print(num_unique_orbits,num_stable_orbits)
    if num_unique_orbits == 0 and num_stable_orbits == 0:
        secondary_counts.append(0)
    else:
        factor = num_stable_orbits/num_unique_orbits
        index = boundaries.index(section)
        secondary_counts.append(count[index]*factor)


xoriginal,yoriginal=plot_piecewise_func(boundaries,count)
plt.plot(xoriginal,yoriginal,'r--',label="$a_{out}(1-e_{out})$ (parent)")
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



plt.plot(xstable,ystable,'r',label="$a_{out}(1-e_{out})$ (stable)")
plt.legend()
plt.xlabel('Periapsis of outer orbit and apoapsis of inner orbit' + r' $(R_{\odot})$')
plt.savefig('Closest Approach distance distributions')




'''For plotting the eccentricity parent and stable distributions'''

df = pd.DataFrame({'eccentricities considered':samples_eccentricities,'eccentricity counts':eccentricity_counter})

count, bins, ignored = plt.hist(samples_eccentricities, 15, normed=True,visible=False)
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
    num_unique_orbits = len(subsection['eccentricities considered'])*len(qoutlist)*len(inclinations)*len(ain)
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
plt.savefig('Eccentricity distributions')





'''For plotting mass distributions'''

df = pd.DataFrame({'mass ratios considered':qoutlist,'mass ratio counts':q_final_counter})

count, bins, ignored = plt.hist(qoutlist, 30, normed=True,visible=False)
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
    num_unique_orbits = len(subsection['mass ratios considered'])*len(ain)*len(samples_eccentricities)*len(inclinations)
    num_stable_orbits = np.sum(subsection['mass ratio counts'])
    if num_unique_orbits == 0 and num_stable_orbits == 0:
        secondary_counts.append(0)
    else:
        factor = num_stable_orbits/num_unique_orbits
        print(factor)
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
plt.savefig('Mass ratio distributions')





'''This is for plotting the original histogram with probability density on the y axis 
and log(a/R(sun)) on the x axis'''
#Convert back to logs for plotting 
ain_final =[]
for i in ain:
    ain_final.append(math.log10(i))

aout_final = []
for i in aout:
    aout_final.append(math.log10(i))
    
    
    
df = pd.DataFrame({'ain_considered':ain_final,'ain_stable_counts':a_ratio_counter})
df2 = pd.DataFrame({'aout_considered':aout_final,'aout_stable_counts':a_ratio_counter})

#AIN
count, bins, ignored = plt.hist(ain_final, 30, normed=True,visible=False)
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
    num_unique_orbits = len(subsection['ain_considered'])*len(qoutlist)*len(samples_eccentricities)*len(inclinations)
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



###AOUT
count, bins, ignored = plt.hist(aout_final, 30, normed=True,visible=False)
bins = list(bins)

boundaries = []
for boundary in bins:
    if boundary == bins[len(bins)-1]:
        pass
    else:
        boundaries.append((boundary,bins[(bins.index(boundary)+1)]))
        
        

secondary_counts = []
for section in boundaries:
    subsection = df2[(df2['aout_considered'] > section[0]) & (df2['aout_considered'] < section[1])]
    #CHECK ASSUMPTION !!!!!!!!!NEED TO CHECK THIS !!!!!--!
    num_unique_orbits = len(subsection['aout_considered'])*len(qoutlist)*len(samples_eccentricities)*len(inclinations)
    num_stable_orbits = np.sum(subsection['aout_stable_counts'])
#    print(num_unique_orbits,num_stable_orbits)
    if num_unique_orbits == 0 and num_stable_orbits == 0:
        secondary_counts.append(0)
    else:
        factor = num_stable_orbits/num_unique_orbits
        index = boundaries.index(section)
        secondary_counts.append(count[index]*factor)
#    print(section,count[index])

xoriginal,yoriginal=plot_piecewise_func(boundaries,count)
plt.plot(xoriginal,yoriginal,'r--',label="$a_{out}$ (parent)")
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




plt.plot(xstable,ystable,'r',label="$a_{out}$ (stable)")
plt.legend()
plt.xlabel(r'$log(\frac{a}{R_\odot}$)')
plt.savefig('Orbital seperation distributions')

                             

                            
                         
                            
                            
                            
                            
                            
                            
                            
                    