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

### Morphology parameters
If you want to use a plasticity model (which for now requires the usage of neuron_type HeteroLIFNeuron or the poisson equivalent, synapse_type HeteroCurrentSynapse and connectivity_type HeteroConnectivity), you will need to specify extra parameters in the neuron section (These can be seen in tests 10 and 12). These parameters determine the morphology/physiology of the neuron population in relation to synaptic weights only.
- morphology_type: Plasticity model that is used. Depending on the model (shown in the [complete guide]), extra parameters will need to be specified.
- dendritic_length/morphology_branch_length: Determines of the theoretical dendrite or length of each branch for each neuron in micrometers. Its ratio with synaptic_gap defines the total spaces or slots available for synapse allocation. Trying to allocate more synapses than available slots will result in a runtime exception.
- synaptic_gap: Determines the space between synapse slots.
Implemented model only work in excitatory populations/synapses. Usage in inhibitory populations will ignore the plasticity framework and the weights will remain constant or the weights will vary but be positive depending on the model unless the opposite is explicitly stated.

### BranchedMorphology parameters
In case the morphology model used works on branched dendrites, the following parameters will apply:
- dendrite_branchings: Number of branchings in the dendritic tree. The generated tree is a binary one, so the number of branches will be equal to $2^{branchings}$.
- synapse_allocation: It can be either 'random' or 'ordered'. It corresponds to the order in which available synapse slots are allocated when the connectivity class connects a presynaptic population. **If connectivity is randomized, the resulting allocation will be randomized no matter what option you select**

## Stimulus Parameters
- type: The type of stimulus which is fed in the network. Here, it is set to WhiteNoiseStimulus. This stimulus is the most standard one, and neurons are fed a white noise input. This input type has two parameters: meanCurrent and sigmaCurrent. Each of these can be defined as a succession of steps in time. 
- meanCurrent: the mean input to each population during a specific step. Here, there are three steps: one step until 4s; another step between 4 and 5s and one after 5s until the end of the simulation. Both population get the same input : 100mV/s during the first and third step and only 10mV/s during the second step.
- sigmaCurrent: the standard deviation of the external input. Here, there is only one step for the whole simulation and the standard deviation is 1mV/s for both populations.

## Recorder Parameters
Each line is a type of output file that can be generated from the simulation. More information is available in [Output Files](README_OutputFiles.md). In the default case, all values are set to 0, so only the default Data file is created. The parameter binSize determines the size of the bins with which time is binned. Here, the value is 10 ms, so the Data file will report averages over 10ms of various states of the system.

## Synaptic Parameters
Every pair of population has its own Synapse object. Here, we have the four: 0<-0, 1<-0, 0<-1 and 1<-1 (first number is postsynaptic population, second presynaptic). Each of them is defined separately.

- type: the model used for the synapse. Here, all four are set to Current Synapse, which is the simplest case: current-based synapses, without plasticity and the PSP is applied at a single time point.
- connected: either 'true' or 'false'. Allows you to skip a certain Synapse between two populations, and avoid its required computations during the simulation.
- D_min and D_max: these parameters define the range of synaptic delays. Each synapse has its own synaptic delay D which is randomly generated at the start of the simulation. D is randomly distributed between D_min and D_max. Here, D_min and D_max are both set to 0 for all synapses, so the PSP is applied at the same time step as the presynaptic spike in all synapses. 
- J: the synaptic strength. This is the value of the PSP (the change in potential on the postsynaptic neuron upon presynaptic spike). Note that this is the only point at which we define which population is excitatory and which is inhibitory. Here, the population 0 is E and the population 1 is I: we set all EPSPs to 0.001 mV and IPSPs to -0.005mV.
- Jpot and Ppot: Some synapses can be potentiated. With a probability of Ppot, synapses have a Synaptic strength of Jpot instead of J. Here, Ppot is set to 0 in all synapses, so no synapse is potentiated.
- connectivity type: the rule for neuron connections. Here, it is set to RandomConnectivity in all synapses. This connectivity rule ensures that each neuron gets the same number of connections, which is determined from the probability of connection and the number of neurons in the presynaptic population: $C=p*N$. For each presynaptic neuron, its presynaptic connections are randomly drawn from the presynaptic population. Here, the connection probability (ConnectProba) is set to 5% in all synapse.
### Branched parameters
In case the morphology model used works on branched dendrites, the following parameters will apply:
- target_branch: It can be either 'ordered' (starts allocation at branch 0), 'random', or a specific branch ID (from zero onwards). In the future, subregion implementation will allow to target multiple branches. 

##Iterative parameters
In every Parameters.txt file you have the option of including the iterate_1 and iterate_2 parameters (recommended to add at the end for consensus). Its usage can be seen in Test10_Parameters.txt, but because of how the programme iterates over these parameters, they will not be visible in the reconstituted copy of the parameter file in the output folder. Instead, it can only be recovered in the input_copy_parameters.txt, which is an identical copy of the original file.
###Syntax example
->iterate_1 neuron_0_noNeurons 50 100 200 400
->iterate_1 synapses_1_0_J 0.05 0.025 0.0125 0.00625
->iterate_2 neuron_1_noNeurons 400 200 100

**Iterate parameters of the same number (1 or 2) must have the same number of entries/columns** 
Not following this recommendation will result in undefined behaviour or crashes. 

Iterate parameters of the same number will be combined by column or number of entry. In this example, 50 neurons will be paired with 0.05 dmV/s, and the next time parameters change, the pair will be 100 neurons with 0.025 dmV/s.

**The original entry in the parameter that is substituted by the iterate_n entries will be ignored**
###Functionality
From the example above, once you run the programme, it will start a simulation with the first column of parameters for the specified parameters. Once it's done, it will start a new simulation with the first column of all iterate_1 parameters and the second column of iterate_2. Once all entries of iterate_2 have been an input in a simulation, the iterate_1 parameters will jump to the second entry, an start anew with the iterate_2 parameter entries. 
The execution will finish once all possible combinations between columns of iterate_1 and iterate_2 have been run. The output writing is independent in each individual simulation, so in the case of crash or abort the comnpleted simulations will have intact data (including total runtime in seconds).

This parameteris extremely effective in case it is necessary to perform a grid search of certain parameters. This implementation is limited to iterate_1 and iterate_2 because of the extreme length resulting from a single run. 

##ParameterOptions
The ParameterOption file is a catalogue of existing classes and can be used to check the other available classes for Stimulus, Neuron, Synapse and connectivity. It also shows the parameters which must be defined for each, and the recognized syntax.

