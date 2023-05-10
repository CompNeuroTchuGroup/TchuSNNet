#ifndef NeuralNetwork_HPP
#define NeuralNetwork_HPP

#include <vector>
#include <fstream>
#include <stdio.h>
#include <string>
#include <iostream>
#include <cstring>
#include <chrono>
#include <random>
#include <ctime>
#include <sys/stat.h>
//#include <unistd.h>
#include "Recorder.hpp"
#include "AdvancedRecorder.hpp"
#include "SynapseSample.hpp"
#include "GlobalFunctions.hpp"
#include "NeuronPopSample.hpp"
#include "Synapse/Synapse.hpp"
#include "Stimulus/UncorrelatedPoissonLikeStimulus.hpp"
#include "Stimulus/WhiteNoiseStimulus.hpp"
#include "Stimulus/WhiteNoiseRescaled.hpp"
#include "Stimulus/WhiteNoiseLinear.hpp"
#include "Stimulus/SpatialGaussianStimulus.hpp"
#include "Stimulus/SpatialPoissonStimulus.hpp"

class NeuralNetwork
{
private:

    GlobalSimInfo  info; //Simulation Parameters

    NeuronPopSample       *neurons;
    SynapseSample         *synapses;
    Recorder              *recorder;
    Stimulus              *stimulus;

    void SaveParameters();
    int  LoadParameters(std::string baseDir,std::vector<FileEntry> *parEntries);
    int  WellDefined();

    void SaveParameterOptions();
public:
    NeuralNetwork(std::string baseDir,std::vector<FileEntry> *parEntries);//std::string baseDir);
    ~NeuralNetwork()=default;

    int  Simulate();

    void makeInputCopy(const std::string&);
    void outputHeteroEvents();
    AdvancedRecorder& GetRecorder(){return *dynamic_cast<AdvancedRecorder*>(recorder);}
};

#endif // NeuralNetwork_HPP
