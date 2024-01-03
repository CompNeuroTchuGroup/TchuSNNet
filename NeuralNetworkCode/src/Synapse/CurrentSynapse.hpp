#ifndef CURRENTSYNAPSE_HEADER_
#define CURRENTSYNAPSE_HEADER_
#include "../GlobalFunctions.hpp"
#include "../NeuronPop/NeuronPop.hpp"
#include "Synapse.hpp"
#include <iostream>
#include <random>
#include <vector>

class CurrentSynapse : public Synapse {
protected:
  std::vector<double> AdvectSpikers(NeuronInt spiker) override;

public:
  CurrentSynapse(PopPtr targetPop, PopPtr sourcePop, GlobalSimInfo *infoGlobal);

  //*****************************
  //******* Get Functions *******
  //*****************************
  int                 GetNoDataColumns() const override { return 1; }
  std::string         GetDataHeader(int dataColumn) override;
  std::string         GetUnhashedDataHeader() const override;
  std::vector<double> GetSynapticState(NeuronInt sourceNeuron) const override;
  std::string         GetTypeStr() const override { return IDstringCurrentSynapse; }

  void SaveParameters(std::ofstream &wParameterStream, std::string idString) const override;
  void LoadParameters(const std::vector<FileEntry> &parameters) override;
};

#endif // CURRENTSYNAPSE
