#ifndef _POWER_LAW_SYNAPSE_HEADER_
#define _POWER_LAW_SYNAPSE_HEADER_

#include <iostream>
#include <random>
#include <vector>

#include "../GlobalFunctions.hpp"
#include "../NeuronPop/NeuronPop.hpp"
#include "Synapse.hpp"

class PowerLawSynapse : public Synapse {
protected:
  double                           kExponent_PLaw{1.0};
  int                              noSpikesForFiringRate{2}; // Number of spikes used to average ISI and invert for average frequency
  std::vector<int>                 spikeCountVector;
  std::vector<std::vector<double>> tableISI;

  std::vector<double> AdvectSpikers(NeuronInt spiker) override;
  // void advect_finalize(std::vector<std::vector<double>> * waiting_matrix) override {};

public:
  PowerLawSynapse(PopPtr targetPop, PopPtr sourcePop, GlobalSimInfo *infoGlobal);
  ~PowerLawSynapse() override = default;

  //*****************************
  //******* Get Functions *******
  //*****************************
  int                 GetNoDataColumns() const override { return 1; }
  std::string         GetDataHeader(int dataColumn) override;
  std::string         GetUnhashedDataHeader() const override;
  std::vector<double> GetSynapticState(NeuronInt sourceNeuron) const override;
  std::string         GetTypeStr() const override { return IDstringPowerLawSynapse; };

  void SaveParameters(std::ofstream &wParameterStream, std::string idString) const override;
  void LoadParameters(const std::vector<FileEntry> &parameters) override;
};

#endif // _POWER_LAW_SYNAPSE_HEADER_
