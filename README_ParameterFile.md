# Parameter File

The parameter file contains all the parameters needed to run a simulation. An example parameter file can be found in the main folder. 

Each entry of the parameter file is a line. The commented lines (#) are not read by the code. The line starts with the name of the parameter, and followed by its value (either numeric, or the name of the referred class). 

As an example, here is what the example parameter file defines:

## General parameters
- Title: the title is "Simulation_Title". This is also going to be the name of the folder containing the output of the simulation.
- SimulationTime: the duration of the simulation is here 10s. Note that the units are only provided for information and cannot be changed: The "seconds" is not read by the code when the parameters are imported.
- dt: The time step of the simulation, here set to 0.01 ms. As a rule of thumb, we observed that using a time step that is 1/1000 of the membrane time constant gives a good tradeoff between accuracy and simulation speed.
- global seed: the seed of the random number generator, used to generate the connectivity and the input noise. Must be an integer.
- Spatial parameters: only relevant if the network has a spatially distributed input or a distance-dependant connectivity (SpatialGaussianStimulus, SpatialPoissonStimulus, DistanceConnectivity). If these classes are used, the number of dimensions must be set under "Dimensions". The size of the network is determined by the spatial density of neurons. In the present example, both are set to 0 since the network does not have any spatial structure.
- Scaling parameters: This option is used to scale the synaptic strength and external stimulus with the network size. The input strength is scaled by N^s, where N is the total number of neurons, and s is the parameter scalingSynapticStrength. In the default case, without scaling, s is set to 0. Typically, s could be set at -0.5, to match the 1/sqrt(N) scaling applied in the Balanced State framework. the other parameter (scaling_C_N) is used to determine whether the input is scaled with N^s (set to 1) or C^s (set to 0), where C is the number of presynaptic connections each neuron has.

## Neuron Parameters
- noPopulations: the number of neuronal populations. It is here set to 2. Each population is homogeneous, and the neuronal properties of each population are defined next. 
- noNeurons: the number of neurons in each population. Population_0 has 3000 neurons and Population_1 has 1000.
- type: The neuron model. Here, both are Leaky Integrate and Fire neurons (LIFNeuron). For LIFneuron, 5 parameters must be defined (tauM, vReset, vThresh ,refractoryTime and resetType). If another neuron type was used, different parameters would need to be provided in the Parameter File.
- tauM: the time constant of the neurons. Here it is set to 10 ms for both populations.
- vReset: The membrane potential at which neurons are reset after they fired. Note that the resting potential is set at 0. Here both populations have the reset potential at 0mV.
- vThresh: the firing threshold, at which neurons fire. Here it is set to 1mV in both populations.
- refractoryTime: the refractory period after firing during which the neuron does not integrate its input and cannot fire. Here it is set to 0s in both populations.
- resetType: this parameter determines whether the excess depolarization after spiking is stored. If set to 0 ,there is a "hard reset" to vReset. If set to 1, the overshoot in membrane potential will be added to the reset potential. Here, both populations are of reset type 0 (hard reset).

## Stimulus Parameters
- type: The type of stimulus which is fed in the network. Here, it is set to WhiteNoiseStimulus. This stimulus is the most standard one, and neurons are fed a white noise input. This input type has two parameters: meanCurrent and sigmaCurrent. Each of these can be defined as a succession of steps in time. 
- meanCurrent: the mean input to each population during a specific step. Here, there are three steps: one step until 4s; another step between 4 and 5s and one after 5s until the end of the simulation. Both population get the same input : 100mV/s during the first and third step and only 10mV/s during the second step.
- sigmaCurrent: the standard deviation of the external input. Here, there is only one step for the whole simulation and the standard deviation is 1mV/s for both populations.

## Recorder Parameters
Each line is a type of output file that can be generated from the simulation. More information is available in the Output Files Wiki. In the default case, all values are set to 0, so only the default Data file is created. The parameter binSize determines the size of the bins with which time is binned. Here, the value is 10 ms, so the Data file will report averages over 10ms of various states of the system.

## Synaptic Parameters
Every pair of population has its own Synapse object. Here, we have the four: 0→0, 1→0, 0→1 and 1→1. Each of them is defined separately.

- type: the model used for the synapse. Here, all four are set to Current Synapse, which is the simplest case: current-based synapses, without plasticity and the PSP is applied at a single time point.
- D_min and D_max: these parameters define the range of synaptic delays. Each synapse has its own synaptic delay D which is randomly generated at the start of the simulation. D is randomly distributed between D_min and D_max. Here, D_min and D_max are both set to 0 for all synapses, so the PSP is applied at the same time step as the presynaptic spike in all synapses. 
- J: the synaptic strength. This is the value of the PSP (the change in potential on the postsynaptic neuron upon presynaptic spike). Note that this is the only point at which we define which population is excitatory and which is inhibitory. Here, the population 0 is E and the population 1 is I: we set all EPSPs to 0.001 mV and IPSPs to -0.005mV.
- Jpot and Ppot: Some synapses can be potentiated. With a probability of Ppot, synapses have a Synaptic strength of Jpot instead of J. Here, Ppot is set to 0 in all synapses, so no synapse is potentiated.
- connectivity type: the rule for neuron connections. Here, it is set to RandomConnectivity in all synapses. This connectivity rule ensures that each neuron gets the same number of connections, which is determined from the probability of connection and the number of neurons in the presynaptic population: C=p*N. For each presynaptic neuron, its presynaptic connections are randomly drawn from the presynaptic population. Here, the connection probability (ConnectProba) is set to 5% in all synapse.


The ParameterOption file is a catalogue of existing classes and can be used to check the other available classes for Stimulus, Neuron, Synapse and connectivity. It also shows the parameters which must be defined for each, and the recognized syntax.