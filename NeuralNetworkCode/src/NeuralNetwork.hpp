#ifndef NEURAL_NETWORK_HPP
#define NEURAL_NETWORK_HPP

#include "SynapseSample.hpp"
#include "NeuronPopSample.hpp"
#include "GlobalFunctions.hpp"
#include "Stimulus/UncorrelatedPoissonLikeStimulus.hpp"
#include "Stimulus/WhiteNoiseStimulus.hpp"
#include "Stimulus/WhiteNoiseLinear.hpp"
#include "Stimulus/SpatialGaussianStimulus.hpp"
#include "Stimulus/SpatialPoissonStimulus.hpp"
#include "Recorder.hpp"
#include <sys/stat.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <string>
#include <iostream>
#include <cstring>
#include <chrono>
#include <random>
#include <ctime>

class NeuralNetwork {
private:

    GlobalSimInfo  infoGlobal; //Simulation Parameters

    std::shared_ptr<NeuronPopSample>       neurons;
    std::shared_ptr<SynapseSample>         synapses;
    std::shared_ptr<Recorder>              recorder;
    std::shared_ptr<Stimulus>              stimulus;

    void SaveParameters();
    void  LoadParameters(std::string baseDirectory,std::vector<FileEntry>& parameterEntries);
    bool  WellDefined();

    void SaveParameterOptions();
public:
    NeuralNetwork(std::string baseDirectory,std::vector<FileEntry> parameterEntries);
    ~NeuralNetwork()=default;

    void  Simulate();

    void makeInputCopies(const std::string&);
    Recorder& GetRecorderReference(){return *recorder;}
};

#endif // NeuralNetwork_HPP
