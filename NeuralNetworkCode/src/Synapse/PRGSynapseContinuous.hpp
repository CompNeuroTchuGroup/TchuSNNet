#ifndef PRGSYNAPSECONTINUOUS_HPP
#define PRGSYNAPSECONTINUOUS_HPP

#include "../GlobalFunctions.hpp"
#include "../NeuronPop/NeuronPop.hpp"
#include "MongilloSynapseContinuous.hpp"
#include "Synapse.hpp"
#include <iostream>
#include <random>
#include <vector>

class PRGSynapseContinuous : public MongilloSynapseContinuous {
protected:
  std::vector<std::vector<double>> l_PRG;

  double tauL_PRG{};      // Decay time constant
  double M_PRG{};         // LPA2 filling probability per transmitted spike
  double deltaU_PRG{};    // Effect of LPA on U
  double deltaTauF_PRG{}; // Effect of LPA on tauF_Mongillo

  std::vector<double> AdvectSpikers(NeuronInt spiker) override;
  double              TransmitSpike(NeuronInt targetId, NeuronInt spikerId) override;
  void                ConnectNeurons() override;

public:
  PRGSynapseContinuous(PopPtr targetPop, PopPtr sourcePop, GlobalSimInfo *infoGlobal);
  ~PRGSynapseContinuous() override = default;

  int                 GetNoDataColumns() const override { return 5; } // J, <y>, <x>, <l>, <submitted_spikes>
  std::string         GetDataHeader(int dataColumn) override;
  std::string         GetUnhashedDataHeader() const override;
  std::vector<double> GetSynapticState(NeuronInt sourceNeuron) const override;
  std::string         GetTypeStr() const override { return IDstringPRGSynapseContinuous; }

  void SaveParameters(std::ofstream &wParameterStream, std::string idString) const override;
  void LoadParameters(const std::vector<FileEntry> &synapseParameters) override;
};

#endif // PRGSYNAPSECONTINUOUS_HPP
