#ifndef SYNAPSESAMPLE_HPP
#define SYNAPSESAMPLE_HPP

#include "GlobalFunctions.hpp"
#include "NeuronPopSample.hpp"
#include "Synapse/CurrentSynapse.hpp"
#include "Synapse/ExponentialCurrentSynapse.hpp"
#include "Synapse/MongilloSynapse.hpp"
#include "Synapse/MongilloSynapseContinuous.hpp"
#include "Synapse/PModelSynapse.hpp"
#include "Synapse/PRGSynapseContinuous.hpp"
#include "Synapse/PowerLawSynapse.hpp"

#include <iostream>
#include <random>
#include <string>
#include <vector>

// class NeuralNetwork;

class SynapseSample {
protected:
  // int                                  generalSynapseSeed;
  // int                                  global_D_max;

  GlobalSimInfo                       *infoGlobal;
  std::shared_ptr<NeuronPopSample>     neurons;
  std::vector<std::vector<SynapsePtr>> synapses; // Indexing is first target neuronPop, second is source neuronPop
  std::vector<Synapse *>               connectedSynapses;
  std::vector<std::vector<bool>>       synapseStates;
  PopInt                               totalNeuronPops;

  std::mutex _syndVMutex;

  void LoadParameters(const std::vector<FileEntry> &synapseParameters);
  void SaveSynapseType(std::string parameterName, std::string type, const std::vector<FileEntry> &synapseParameters);
  void SetUpSynapse(PopInt targetPop, PopInt sourcePop, std::string type, std::vector<FileEntry> synapseParameters);

public:
  SynapseSample(std::shared_ptr<NeuronPopSample> neurons, std::vector<FileEntry> &synapseParameters, GlobalSimInfo *infoGlobal);
  ~SynapseSample() = default;

  void ConnectNeurons();
  void WriteConnectivity(std::string filename, NeuronInt noNeuronsConnectivity) const;
  void WriteDistributionD(std::string filename, NeuronInt noNeuronsDelay) const;
  void WriteDistributionJ(std::string filename, NeuronInt noNeuronsJPot) const;

  //*******************
  // Get-Functions
  //*******************
  // int         GetNoDataColumns() const;
  // CAREFUL CALLING THESE FUNCTIONS, ALWAYS CHECK CONNECTED STATE
  int         GetNoDataColumns(PopInt targetPop, PopInt sourcePop) const;
  std::string GetDataHeader(int dataColumn) const;
  std::string GetUnhashedDataHeader() const;
  // int         GetSeed();
  // int         GetMaxD(){return global_D_max;}
  /**
   * Returns the sum of the synaptic state variables of neuron pre_neuron
   * targeting neurons in post_population.
   */
  // CAREFUL CALLING THESE FUNCTIONS, ALWAYS CHECK CONNECTED STATE
  std::vector<double> GetSynapticState(PopInt targetPop, PopInt sourcePop, NeuronInt sourceNeuron) const;
  double              GetRecurrentInput(PopInt targetPop, PopInt sourcePop, NeuronInt targetNeuron) const;
  double              GetCumulatedDV(PopInt targetPop, PopInt sourcePop) const;
  NeuronInt           GetNoTargetedNeurons(PopInt targetPop, PopInt sourcePop, NeuronInt sourceNeuron) const;
  bool                GetConnectedState(PopInt targetPop, PopInt sourcePop) const { return synapseStates.at(targetPop).at(sourcePop); };

  //**************************************
  // Functions which are redirected to the synapse types
  //**************************************
  void Advect(std::vector<std::vector<double>> &synaptic_dV);
  void Reset();
  void SaveParameters(std::ofstream &wParameterStream) const;

  // void Test();
};

#endif // SYNAPSESAMPLE_HPP
