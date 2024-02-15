
#ifndef _NEURONPOP_BASE_CLASS_HPP
#define _NEURONPOP_BASE_CLASS_HPP
// Include order is important to avoid circular dependency problems
#include "./Morphology/HybridModels/AlphaResourceHSTDP.hpp"
#include "./Morphology/HybridModels/HeteroGraupnerBrunel.hpp"
#include "./Morphology/HybridModels/MACRbPModel.hpp"
#include "./Morphology/TimingDependentModels/MonoDendriteSTDPBiWindow.hpp"
#include "./Morphology/TimingDependentModels/MonoDendriteSTDPTazerart.hpp"
#include "./Morphology/TimingDependentModels/MonoDendriteSTDPTazerartRelative.hpp"
#include <execution>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

class NeuronPop {
protected:
  GlobalSimInfo *infoGlobal;

  PopInt    identifier{};
  NeuronInt noNeurons{};
  TStepInt  refractorySteps{};
  double    membraneVTau{};
  double    membraneExpDecay{0.002};
  double    resetV{};
  double    thresholdV{1.0};
  // double targetRate;
  double prevMinFrequency{0.333};
  bool   definedPrevMeanFreq{false};

  int          seed{};
  std::mt19937 generator;
  bool         userSeeds{false};

  bool                                     hasPlasticity{false}, isBranched{false};
  bool                                     savedModel{false};
  std::vector<std::unique_ptr<Morphology>> morphology{};

  std::vector<double>    membraneV;             // membrane potential
  std::vector<NeuronInt> spikerNeurons;         // indices of all neurons that have emitted a spike in the previous time step
  std::vector<NeuronInt> spikerNeuronsPrevdt;   // indices of all neurons that have emitted a spike in the previous time step
  std::vector<TStepInt>  previousSpikeDistance; // last spike time for every neuron
  std::vector<double>    xAxisPosition;         // position on the x axis//No reason for valarray
  std::vector<double>    yAxisPosition;         // position on the x axis
  // std::valarray<long>     global_id;  // global neuron id (unique across populations)
  void ClearSpikerVector();

  std::string morphologyType{"None"};
  // Necessary mutex locks

  std::mutex _connectMutex;

public:
  // constructor
  NeuronPop(GlobalSimInfo *infoGlobal, PopInt popID);
  virtual ~NeuronPop() = default;

  //*******************
  // Get-Functions
  //*******************
  NeuronInt GetNoNeurons() const { return noNeurons; }
  // long GetGlobalNeuronId(long i){return global_id[i];}
  double GetXPosition(NeuronInt neuron) const { return xAxisPosition.at(neuron); }
  double GetYPosition(NeuronInt neuron) const { return yAxisPosition.at(neuron); }

  // those functions return parameterValues and addresses ( pay attention here!!!)
  // TODO/Suggestion: Replace pointer return with a return by constant reference.
  // https://github.com/saiftyfirst/BP_Demos/blob/master/C%2B%2B/constRef_vs_pointer.cpp
  const std::vector<NeuronInt> &GetSpikers() const { return spikerNeurons; }
  const std::vector<NeuronInt> &GetSpikersPrevDt() const { return spikerNeuronsPrevdt; }

  // long                get_previous_spike_step(long i) {return previousSpikeStep[i];}
  double GetTimeSinceLastSpike(NeuronInt neuron) const { return static_cast<double>(previousSpikeDistance.at(neuron)) * infoGlobal->dtTimestep; }
  double GetPotential(NeuronInt neuron) const { return membraneV.at(neuron); }

  virtual std::string GetType() const = 0;
  std::string         GetMorphologyType() const { return morphologyType; }
  PopInt              GetId() const { return this->identifier; }

  SynInt GetNoSynapses() const; // This is purely for virtualization reasons
  int    GetNoBranches() const { return static_cast<BranchedMorphology *>(morphology.at(0).get())->GetNoBranches(); }
  int    GetMaxGapDelay(int delayPerGap) const { return (morphology.at(0).get())->GetMaxGapDelay(delayPerGap); }
  double GetSynapticDistanceToSoma(int neuronId, int synapseId) const { return morphology.at(neuronId)->GetSynapticDistanceToSoma(synapseId); }
  //*******************
  // Set-Functions
  //*******************
  void SetNeurons();
  void SetPosition();
  // void SetGlobalNeuronId(long global_neuron_id);

  // Main methods
  virtual void Advect(const std::vector<double> &synaptic_dV) = 0;
  void         AdvectPlasticityModel();
  // void    LoadParameters(const std::vector<FileEntry>& neuronParameters,signed long noNeurons);
  virtual void LoadParameters(const std::vector<FileEntry> &neuronParameters);
  void         LoadPlasticityModel(const std::vector<FileEntry> &morphologyParameters);
  virtual void SaveParameters(std::ofstream &wParameterStream) const;
  void         SavePlasticityModel(std::ofstream &wParameterStream, std::string idString) const;
  void         ModelSaved() { savedModel = true; }

  // Dyn_casting optimization methods
  bool HasPlasticityModel() const { return hasPlasticity; }
  bool IsBranchedBool() const { return isBranched; }
  bool ignoreJDParameters() const;
  bool HasSteadyState() const;

  // All of the following functions throw to ease the virtualization, but it is a bad coding practice

  void RecordExcitatorySynapticSpike(NeuronInt neuronId, BaseSpinePtr spinePtr) { this->morphology.at(neuronId)->RecordExcitatoryPreSpike(spinePtr); }
  void RecordPostSpike(NeuronInt spikerID) { this->morphology.at(spikerID)->RecordPostSpike(); }

  BaseSpinePtr        AllocateNewSynapse(NeuronInt neuronId, BranchTargeting &targeting);
  std::string         GetIndividualSynapticProfileHeaderInfo() const;
  std::string         GetOverallSynapticProfileHeaderInfo() const;
  std::vector<double> GetIndividualSynapticProfile(NeuronInt neuronId, NeuronInt spineID) const;
  std::vector<double> GetOverallSynapticProfile(NeuronInt neuronId) const;

  std::vector<std::pair<std::string, double>> GetSteadyStateData() const;

  void PostConnectSetUp(); // This function cannot be called without throwing first};
};

#endif // NeuronPop_HPP
