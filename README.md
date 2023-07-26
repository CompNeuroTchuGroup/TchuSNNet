# Welcome to our Spiking network Simulation code "BrainCode"

## Theory of Neural Dynamics Group
### Tatjana Tchumatchenko

Contributors:

- Ahmed, Saif 
- Becker, Sophia
- Bertolin, Antoni
- Dutine, Alexander
- Ekelmans, Pierre
- Konrad, Sara
- Kovalevskaya, Elizaveta (Group of Algorithm Engineering)
- Nold, Andreas 
- Penschuck, Manuel (Group of Algorithm Engineering)

- Theory of neural dynamics, Max-Planck-Institut for Brain Research, Frankfurt
- Theory of neural dynamics, Institute of experimental epileptology and cognition research, University Hospital Bonn
- Theory of neural dynamics, Institute of Physiology, Johannes Gutenberg University Mainz
- Group of Algorithm Engineering, Institut für Informatik, Goethe University Frankfurt

February 19th, 2021


## Wiki Manual
Please check the Wiki manual page on https://gitlab.mpcdf.mpg.de/mpibrano/neuronal-network-simulations/-/wikis/home for more information about the practical aspects of compiling and running the code is accessible in the Overview GitLab Wiki.

## General Information
This code is written in C++, and all its files are contained in the folder NeuralNetworkCode.

All the parameters of a simulation are contained in a parameter file. The parameter file is a txt file with a specific format that will be read by the simulation code. The code can get up to two arguments:

- If two arguments are provided (./main $BASE $InputFile), the first argument is the path where the output of the simulation will be saved and the second argument is the path to the parameter file
- If only one argument is provided (./main $BASE), it is the path to the folder containing the parameter file; and the output of the simulation will be saved there too
- If no argument is provided (./main ), the parameter file must be located in the same folder as the code

## Tests

We provide eight tests to reproduce published results with the code. Please see [here](README_Tests.md) for more details.

## Parameter Files

An example parameter file is provided in the main folder (Parameters.txt).
Each entry of the parameter file is a line. The commented lines (#) are not read by the code.
The line starts with the name of the parameter, and followed by its value (either numeric, or the name of the referred class). The parameter file is broken down into several sections:

- General parameters (title of the simulation, duration, time step,...)
- Number of populations of neurons
- Properties of each population (number of neurons, neuron type, time constant,...)
- Stimulus (type of stimulus, amplitude, duration,...). A stimulus is defined by a succession of steps, which are written over multiple lines.
- Recording options (recording bin size, choice of output files)
- Synapses (synapse type, connectivity type, synaptic strength,...). For N populations, there will be N^2 synapse types.

More information about the structure of the parameter files is accessible in the [Parameter File Readme](README_ParameterFile.md).

The ParameterOption file is a catalogue of existing classes for Stimulus, Neuron types, Synapse types and connectivity. It also provides the syntax recognized by the simulation code.

## Simulation and Ouputs
When the simulation starts, it creates a folder containing the output of the code (which gets written as the simulation runs). The location of this folder depends on the arguments given and it is named after the title of the simulation. It contains a copy of the parameter file, and all the output files.
In the minimum case, there is one output file, named *Title*_Data, which stores the average potential, firing rates and inputs. More output files can be created (raster, connectivity maps, I and V profile,...) if the corresponding options are chosen in the parameter file under "Recording parameter". More information about the output files is accessible in the Output Files wiki.

To get familiar with the Spiking Network Simulation code, a few tests are provided (Test folder). Each of these tests is used to make sure the code runs properly and returns the right output. In every case, the parameter file needed to run the test is provided, as well as a markdownfile showing the expected result of the simulation.

See [here](README_OutputFiles.md) for more detailed information.

## Model explanations
To see the explanation of the implementation of some of the available classes, go to the [Models](Models.md) file.
