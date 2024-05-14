//
// Created by Saif Ahmed on 15.07.21.
//
//
// Refactored by Antoni Bertolin on 14.06.23
//
#ifndef _MONODENDRITE_STDP_TAZERART_HPP
#define _MONODENDRITE_STDP_TAZERART_HPP

#include "MonoDendriteSTDP.hpp"
#include <string>

class MonoDendriteSTDPTazerart : public MonoDendriteSTDP {

protected:
  void UpdateLTP(signed long synId) override;
  void UpdateLTD(signed long synId) override;

  double gLTP(double deltaT) const override;
  double gLTD(double deltaT) const override;

  double muLTP{};
  double sigmaLTP{};
  double alpha{};

  double muLTD{};
  double sigmaLTD{};
  double beta{};

public:
  explicit MonoDendriteSTDPTazerart(GlobalSimInfo *infoGlobal, const std::vector<FileEntry> &morphologyParameters);
  ~MonoDendriteSTDPTazerart() override = default;
  virtual std::string GetType() const override { return IDstringMonoDendriteSTDPTazerart; }

  void SaveParameters(std::ofstream &wParameterStream, std::string neuronIdentificator) const override;
  void LoadParameters(const std::vector<FileEntry> &parameters) override;
  void CheckParameters(const std::vector<FileEntry> &parameters) override;
};

#endif // NEURALNETWORK_MONODENDRITESTDPTAZERART_H
