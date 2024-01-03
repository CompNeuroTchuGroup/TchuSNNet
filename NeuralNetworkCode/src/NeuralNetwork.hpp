#ifndef NEURAL_NETWORK_HPP
#define NEURAL_NETWORK_HPP

#include "GlobalFunctions.hpp"
#include "NeuronPopSample.hpp"
#include "Recorder.hpp"
#include "Stimulus/SpatialGaussianStimulus.hpp"
#include "Stimulus/SpatialPoissonStimulus.hpp"
#include "Stimulus/UncorrelatedPoissonLikeStimulus.hpp"
#include "Stimulus/WhiteNoiseLinear.hpp"
#include "Stimulus/WhiteNoiseStimulus.hpp"
#include "SynapseSample.hpp"
#include <chrono>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <vector>

class NeuralNetwork {
private:
  GlobalSimInfo infoGlobal; // Simulation Parameters

  std::shared_ptr<NeuronPopSample> neurons;
  std::shared_ptr<SynapseSample>   synapses;
  std::shared_ptr<Recorder>        recorder;
  std::shared_ptr<Stimulus>        stimulus;

  void SaveParameters();
  void LoadParameters(std::string baseDirectory, std::vector<FileEntry> &parameterEntries);
  bool WellDefined();

  void SaveParameterOptions();

public:
  NeuralNetwork(std::string baseDirectory, std::vector<FileEntry> parameterEntries);
  ~NeuralNetwork() = default;

  void Simulate();

  void      MakeInputCopies(const std::string &);
  Recorder &GetRecorderReference() { return *recorder; }
};

#endif // NeuralNetwork_HPP
