# HowTo

Here is an overview on how to get the simulation running.

## Step 1: Download 
The whole projected can be downloaded on your local machine as a compressed folder. Alternatively, you can also clone the repository.

## Step 2: Compile 
Once you have the code on your machine, you will need to compile the c++ code to get an executable file.
There are many different ways to do this. Here we give one possible way for each operating system.
### On Mac

In ```Makefile```, change in the second line the compiler to ```CC = clang++``` and if necessary remove the option ```-fopenmp```. Then type "make" in a command line. Prerequisite is that Xcode is installed. 

### On Windows
Using the Visual Studio IDE, start a new project and include all the .hpp and .cpp files into it. Then make sure to be in "Release" mode, and use "Build" to compile the code. Once the code is compiled, go in the folder where the project was created and look into the "Release" folder for an executable file. This is the compiled simulator.
Newer versions of MS-compiler can fail to compile because of safety warnings. To alleviate this issue, you can add a new line "_CRT_SECURE_NO_WARNINGS" into Preprocessor Definitions; under Projects>Properties>C/C++>Preprocessor>Preprocessor Definitions.

Alternatively, the compilation can also be done without using an IDE, by downloading a gcc compiler and calling it in command line.

### On Linux
Along with the hpp and cpp files, the NeuralNetworkCode folder contains the Makefile to compile.
From the location of the code, type "make" in a command line.
```
>make 
```
## Step 3: Set parameters

All the parameters of a simulation are contained in a parameter file. The parameter file is a txt file with a specific format that will be read by the simulation code. For a sample parameter file and details, see the [Parameter file wiki](Documents/Parameter file).
In case of usage of AdjacencyMatrixConnectivity or DictatNeuronPop it will be required to include their required files (Connectivity matrixes or spiker files) in the same folder. Examples of this are in Test1 and Test11.

## Step 4: Run simulation
Once the code has been compiled and the parameter file is ready, the code is ready to be run. 

### On Windows
The code can be run in command line. To do so, open a terminal in the location of the executable file. The simulation can then be launched by calling the executable with two arguments, only one, or no argument:\
- **With two arguments**\
The first argument is the path where the output folder will be created and the second argument is the full path of the parameter file.
```
>Simulation.exe C:\ResultsFolder\ C:\ParameterFolder\InputFile.txt 
```
- **With one argument**\
The argument is the path of the folder containing the parameter file, which must have the name "Parameters.txt". The output of the simulation will be stored in a subfolder there.
```
>Simulation.exe C:\SimulationsFolder
```
- **With no argument**\
The executable file must be in the same folder as the parameter file, named "Parameters.txt". The output of the simulation will be stored in a subfolder there.
```
>Simulation.exe
```

It is also possible to run the simulation from the Visual Studio IDE. The Working directory or command arguments can be set in Properties>Configuration Properties>Debugging. 

Due to the heavy computation needed for large network or long simulations, it is recommended to run the code from a dedicated compute machine for networks larger than 5000 neurons.

### On Linux and Mac
To run the code, enter "./main" followed with the arguments (if any):

- If two arguments are provided (./main $BASE $InputFile), the first argument is path where the output of the simulation will be saved and the second argument is the path to the parameter file 
```
./main /path/to/datafolder/  /path/to/parameterfile 
```
- If only one argument is provided (./main $BASE), it is the path to the folder containing the parameter file, which must be named "Parameters.txt".The output of the simulation will be saved in the same folder.

```
./main /path/to/folder 
```
- You can run the code without argument (./main ), if parameter file, named "Parameters.txt" is in the current directory.

```
./main
```

## Step 5: Run & Analyse Tests 

The folder containing the output of the simulation is is named after the simulation title and the time of the code execution. All the output files are plain txt files. More information about the output files can be found in the [Output Wiki](Documents/Output-Files)


To get familiar with the Spiking Network Simulation code, a few tests are provided in the Test folder (see [Tests](Tests) for details). Each of these tests is used to make sure the code runs properly and returns the right output. In every case, the parameter file needed to run the test is provided, as well as a markdownfile showing the expected result of the simulation.
