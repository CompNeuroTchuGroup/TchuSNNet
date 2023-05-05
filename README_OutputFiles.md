# Output Files

As the simulation runs, data will be written in multiple Output files. These Files are contained in the output folder, according to the argument that was provided upon execution of the simulation. In the minimal case, there will only be one output file: the Datafile. If more Recording options were selected in the parameter file, more files can be created to record more precisely specific aspects of the system.

**The Data File**

The Data file contains different averages of the system over time. The averages are computed over a time period defined as the binSize (provided in the parameter file). These values.
In this file, every line is a time bin and every column is the population average of some quantity : average membrane potential, average firing rate, average incoming external and recurrent current and average total incoming current. The final columns of the file are the quenched fluctuations (variability of the total input within the population) and the temporal fluctuations (variability of the population averaged total input within the time bin).
The header of the file comes with a definition of what each column is.

**Connectivity Matrices**

If the parameter noNeuronsConnectivity is set to a non-zero integer, Connectivity_Matrix files will be created.
These files record the connectivity of the first x neurons (where x is the parameter provided in the parameter file). Given N populations of neurons, N^2 files will be created, corresponding to all pairs of post- and pre- synaptic populations.
In the connectivity Matrix, each row is a presynaptic neuron and each column is a post-synaptic neuron. The matrix holds binary values: 0 if the neuron pair is not connected and 1 if they are.

**Delay Matrices**

Similarly to the Connectivity matrices, the parameter noNeuronsDelay determines whether DelayConnectivity_Matrix files are created, and how many neurons they record. They follow the same structure as Connectivity_Matrix, but instead of providing binary values, they provide the synaptic delay of each synaptic connection, given as a number of time steps. If the neurons are not connected, the value is nan instead.

**Potentiation Matrices**

Just like for the other two types of matrices, the noNeuronsJPot parameter determines whether JPotConnectivity_Matrix files are created, and how many neurons they record. These matrices show the synaptic strength associated with each of the recorded synapses, 
given in mV/spike. If the neuron pair is not connected, the value is nan instead.

**RasterPlot**

The parameter noRasterPlotNeurons has multiple entries: one entry corresponding to each population and determining how many neurons are recorded from in each population. The last entry is the time at which the recording starts.
If at least one neuron is recorded from, a Rasterplot file will be created, which keeps track of the spike time of the recorded neurons. In this file, the first column is the spike time, the second column is the id of the spiking neuron and the third column is the population the neuron belongs to. 
Note that the neurons which are recorded from are picked across the whole population. The neuron_id are therefore uniformly distributed in the whole range from 0 to the population size.

**NeuronProfiles**

The parameter notrackNeuronProfiles has N entries (one per population) and determines the number of neurons from each population which are going to get recorded. If any neuron is recorded from, two files are created: Potential and Current. The Potential file records the membrane potential (in mV) of neurons at every time step; the Current file records the input that neurons receive (in mV/s). The neurons which are recorded are the first x_i neuron of each population, where x_i the _i_th parameter entry of notrackNeuronProfiles.

**CurrentContribution**

The parameter CurrentContributions works like the noRasterPlotNeurons parameter, with N+1 entries being the number of neurons recorded from in each population and the time at which recording starts. Here again, the recorded neurons are picked across the whole population.
If any neuron is recorded from, a file CurrentContributions is created which records the input to neurons, separated by source. The input to neurons is broken into external feedforward input (IFFd) and Recurrent input (Irec) and the recurrent inputs coming from different populations are separate columns in the file. 
Here again, the header clarifies what each column is.

**Synapse recorder**

The parameter trackSynapses is binary; if it is set to 1, the synapses are recorded from and a file Synapse is created. For each Synapse group, between all pairs of populations, the average state of the synapses is recorded. The specific content of the file depends on the synapse types used, and is mostly used in case of plastic synapses. In the minimal case, the average synaptic strength of the synapses is recorded in each time bin.
The content of the file is clarified by its header.

**Heatmap**

The parameter Heatmap is used for networks that have a spatial structure. If it is set to a non zero integer, a HeatmapRate file will be created for every population. 
If the simulation is run in 1 dimension, the linear dimension of the system is binned into n spatial bins, where n is the parameter Heatmap (provided in the parameter file). In each time bin, the average firing rate in each spatial bin is written in the file.
If the simulation is run in 2 dimensions, the system is split in n^2 spatial bins, and each line of the Heatmap file contains the n^2 firing rates.

**NeuronProfiles**

The parameter notrackNeuronProfiles has N entries (one per population) and determines the number of neurons from each population which are going to get recorded. If any neuron is recorded from, two files are created: Potential and Current. The Potential file records the membrane potential (in mV) of neurons at every time step; the Current file records the input that neurons receive (in mV/s). The neurons which are recorded are the first x_i neuron of each population, where x_i the _i_th parameter entry of notrackNeuronProfiles.

**HeteroSynapses**

The parameter notrackHeteroSynapseProfiles has N+1 entries where N is the number of populations. Each entry will determine how many synapses from each of the neurons tracked in NeuronProfiles are tracked. The file is structured as following: one time column and one column per synapse, with its ID on the top. Depending on the model, the output will differ for each datapoint, which will normally be structured as a tuple of numbers. The contents of this tuple will be written in one of the commented lines of the file at the top.

**OverallHS**

For every neuron that has at least one tracked synapse, the average weight, number of plasticity events and total presynaptic spikes and postsynaptic spikes will be inside a tuple for each datapoint.

**ParsedDataFiles**

For some files, their contents might be too large to store or work with, so when turning parsing ON vs OFF, at the end of the simulation the datafile parser will create a file (for every implemented parsing) that contains the same information but more compact. These files will also contain important metadata, like the number of neurons per population, mangitude of dt in the simulation, number of total timesteps run, population IDs and the simulated time in seconds. To manipulate these files there are special python notebooks in the ./Test_results/AnalyseCode/ folder.

