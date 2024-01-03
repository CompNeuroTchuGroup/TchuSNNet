#ifndef MONGILLOSYNAPSECONTINUOUS_HPP
#define MONGILLOSYNAPSECONTINUOUS_HPP
#include "../GlobalFunctions.hpp"
#include "../NeuronPop/NeuronPop.hpp"
#include "Synapse.hpp"
#include <iostream>
#include <random>
#include <vector>

class MongilloSynapseContinuous : public Synapse {
protected:
  double u_MongilloC{};
  double tauF_MongilloC{}; // in units of time steps
  double tauD_MongilloC{}; // in units of time steps

  std::vector<std::vector<double>> x_MongilloC, y_MongilloC, spikeSubmitted;

  // int seed;
  // std::mt19937 generator;
  // std::uniform_real_distribution<double> uniformDistribution;

  std::vector<double> AdvectSpikers(NeuronInt spiker) override;

  virtual double TransmitSpike(NeuronInt targetId, NeuronInt spikerId);

public:
  MongilloSynapseContinuous(PopPtr targetPop, PopPtr sourcePop, GlobalSimInfo *infoGlobal);
  ~MongilloSynapseContinuous() override = default;

  void ConnectNeurons() override;

  //*****************************
  //******* Get Functions *******
  //*****************************
  int                 GetNoDataColumns() const override { return 4; } // J, <y>, <x_Mongillo>, <submitted_spikes>
  std::string         GetDataHeader(int dataColumn) override;
  std::string         GetUnhashedDataHeader() const override;
  std::vector<double> GetSynapticState(NeuronInt souceNeuron) const override;
  std::string         GetTypeStr() const override { return IDstringMongilloSynapseContinuous; }

  void SaveParameters(std::ofstream &wParameterStream, std::string idString) const override;
  void LoadParameters(const std::vector<FileEntry> &synapseParameters) override;
};

#endif // MONGILLOSYNAPSECONTINUOUS_HPP
