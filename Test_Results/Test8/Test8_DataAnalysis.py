# # Test8
# #### Pierre Ekelmans (translated by Antoni Bertolin)

# This version of the pipeline is unable to run because the convolve2d from scipy is exceedingly slow. Try using the MATLAB version

# ## General Introduction
# This test is designed to reproduce the figure 3b and 4b from the Rosenbaum 2017 paper: 
#     Rosenbaum, Robert, et al. "The spatial structure of correlated neuronal variability." Nature neuroscience 20.1 (2017): 107-114.
# These figures illustrate the covariance between different types of input current depends on distance. In particular, it shows EI balance breaks as the recurrent connections extend further than feedforward connections.
# This test therefore relies on spatially structured networks and the measurement of different types of input current. The features tested are : Distance Connectivity, SpatialPoissonStimulus, EIFNeuron and CurrentContribution recording.


import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import convolve2d




# Import Data
base_path = os.getcwd()
project_name1 = 'Test8A'  # For figure 3B
project_name2 = 'Test8B'  # For figure 4B

# Import data for Fig3B
folder_list = [folder_name for folder_name in os.listdir(base_path) if folder_name.startswith(project_name1)]
folder_name = folder_list[-1]
#os.chdir(os.path.join(base_path, folder_name))

for file in os.listdir(folder_name):
    if file.endswith("CurrentContribution.dat"):
        filename=os.path.join(folder_name, file)
print(filename)

delimiter = '\t'
start_row = 11
format_el = '%f'
format_spec = format_el  # The first column is time
for _ in range(500):
    format_spec += format_el + format_el + format_el  # Each recorded neuron has 3 sources of inputs: recurrent E and I and feedforward
format_spec += '%[^\n\r]'

with open(filename, 'r') as file:
    file_data = file.readlines()
    file.close()

dataArray = []
for line in file_data[start_row-1:]:
    line_data = line.strip().split(delimiter)
    formatted_data = [float(data) for data in line_data]
    dataArray.append(formatted_data)

Test8A = np.array(dataArray)
filename = None  # Clear the filename variable
del dataArray
# Import data for Fig4B
folder_list = [folder_name for folder_name in os.listdir(base_path) if folder_name.startswith(project_name2)]
folder_name = folder_list[-1]
#os.chdir(os.path.join(base_path, folder_name))

for filena in os.listdir(folder_name):
    if filena.endswith("CurrentContribution.dat"):
        filename=os.path.join(folder_name, filena)
print(filename)

with open(filename, 'r') as file:
    file_data = file.readlines()
    file.close()

dataArray = []
for line in file_data[start_row-1:]:
    line_data = line.strip().split(delimiter)
    formatted_data = [float(data) for data in line_data]
    dataArray.append(formatted_data)

Test8B = np.array(dataArray)
filename = None  # Clear the filename variable
del dataArray



Ntot=4e4
Nrec=500
Nx=np.sqrt(Ntot)
Index=np.arange(0, Ntot, Ntot/Nrec)
Irecord=np.array([1/Nx*np.floor(Index/Nx), 1/Nx*np.mod(Index,Nx)]) #(xy coordinates of each neuron recorded)

dt=0.5

# Function to compute distances
distfun=lambda x1,y1,x2,y2 : np.sqrt(np.minimum(np.abs(x1-x2),1-np.abs(x1-x2))**2 + np.minimum(np.abs(y1-y2),1-np.abs(y1-y2))**2)

# distance bins
binsize=0.025
bd=np.append(1/Nx-1e-9, np.arange(binsize, .708, binsize))

# Time constant of filter
tauKc=15
# Low-pass filter
Kc=np.exp(-np.abs(np.arange(-5*tauKc, dt + 5*tauKc, dt))/tauKc)
Kc=Kc/np.sum(Kc)

for it in [0, 1]: #Respectively the 3B and 4B figure; both share the same pipeline
    if it==0:
        Irec=Test8A[1:,1:501].T+Test8A[1:,501:1001].T
        Iffd=Test8A[1:,1001:1501].T
        Tburn=np.round(2/(Test8A[1,0]-Test8A[0,0]))
    else:
        Irec=Test8B[1:,1:501].T+Test8B[1:,501:1001].T
        Iffd=Test8B[1:,1001:1501].T
        Tburn=np.round(2/(Test8B[1,0]-Test8B[0,0]))
    # Two-dim filter


Kcc=np.zeros([np.size(Iffd,0)*2+1,np.size(Kc)])
Kcc[np.size(Iffd,0),:]=Kc
# Low-pass filter currents
Irec0=convolve2d(Irec,Kcc,'same')
del Irec
print("Irec convolved")
Iffd0=convolve2d(Iffd,Kcc,'same')
del Kcc, Iffd
print("Iffd convolved")
    # Get rid of beginning and end, which are corrupted by initial transient and by filtering

Irec0=Irec0[:,int(Tburn):]
Iffd0=Iffd0[:,int(Tburn):]
# locations of recorded neurons
xlocs=Irecord[0,:]
ylocs=Irecord[1,:]
# Distances between neurons
x1,x2=np.meshgrid(xlocs,xlocs)
y1,y2=np.meshgrid(ylocs,ylocs)
print("meshgrids done")
distances=distfun(x1.ravel(),y1.ravel(),x2.ravel(),y2.ravel())
print("distances calculated")
# Compute covariance matrix between all ffwd and all rec inputs


AllCovs=np.cov(np.concatenate((Iffd0, Irec0), axis=0))
print("AllCovs done")
del Irec0, Iffd0
# Get ffwd-rec covariances
FRCovs=AllCovs[0:Nrec,Nrec:]
FRCovs=FRCovs.ravel()
print("FRCovs done")
# Get ffwd-ffwd covs
FFCovs=AllCovs[0:Nrec,0:Nrec]
FFCovs=FFCovs.ravel()
print("FFCovs done")
# Get rec-rec covs
RRCovs=AllCovs[Nrec:,Nrec:]
RRCovs=RRCovs.ravel()
print("RRCovs done")
del AllCovs


TotalCovs=FFCovs+RRCovs+2*FRCovs
print("TotalCovs done")
# Compute mean and stderr of all covariances over each distance bin
I = np.digitize(distances,bd)
mFRCovs, errFRCovs, mFFCovs, errFFCovs, mRRCovs, errRRCovs, mTotalCovs, errTotalCovs = [np.zeros(len(bd)-1) for _ in range(8)]
print("Pre-loop done")
for j in range(len(bd)-1):
    mFRCovs[j]=np.mean(FRCovs[I==j])
    errFRCovs[j]=np.std(FRCovs[I==j])/np.sqrt(np.count_nonzero((I==j)))
    mFFCovs[j]=np.mean(FFCovs[I==j])
    errFFCovs[j]=np.std(FFCovs[I==j])/np.sqrt(np.count_nonzero((I==j)))
    mRRCovs[j]=np.mean(RRCovs[I==j])
    errRRCovs[j]=np.std(RRCovs[I==j])/np.sqrt(np.count_nonzero((I==j)))
    mTotalCovs[j]=np.mean(TotalCovs[I==j])
    errTotalCovs[j]=np.std(TotalCovs[I==j])/np.sqrt(np.count_nonzero((I==j)))
print("Loop done")

bd=bd[:-1]
plt.figure()
plt.plot(bd,mFRCovs/mFFCovs[0],'m',linewidth=2)
plt.plot(bd,mFFCovs/mFFCovs[0],'b',linewidth=2)
plt.plot(bd,mRRCovs/mFFCovs[0],'r',linewidth=2)
plt.plot(bd,mTotalCovs/mFFCovs[0],'k',linewidth=2)
plt.xlabel('distance (a.u.)')
plt.ylabel('Current covariance (normalized)')
plt.legend(['ffd-rec','ffd-ffd','rec-rec','Tot'])
plt.tick_params(labelsize=16)
plt.title('3B' if it==0 else '4B')
plt.axis('tight')
plt.show()