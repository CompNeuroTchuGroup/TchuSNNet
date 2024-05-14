//
// Created by Saif Ahmed on 15.07.21.
//
//
// Refactored by Antoni Bertolin on 14.06.23
//
#ifndef _MONODENDRITE_SYMMETRIC_STDP_HPP
#define _MONODENDRITE_SYMMETRIC_STDP_HPP

#include "MonoDendriteSTDP.hpp"
#include <string>

class MonoDendriteSTDPBiWindow : public MonoDendriteSTDP {

protected:
  void UpdateLTP(signed long synId) override;
  void UpdateLTD(signed long synId) override;

  double gLTP(double deltaT) const override;
  double gLTD(double deltaT) const override;

  double tauLTP{};

  double tauLTD{};

public:
  explicit MonoDendriteSTDPBiWindow(GlobalSimInfo *infoGlobal, const std::vector<FileEntry> &morphologyParameters);
  ~MonoDendriteSTDPBiWindow() override = default;
  virtual std::string GetType() const override { return IDstringMonoDendriteSTDPBiWindow; }

  void SaveParameters(std::ofstream &wParameterStream, std::string neuronIdentificator) const override;
  void LoadParameters(const std::vector<FileEntry> &parameters) override;
  void CheckParameters(const std::vector<FileEntry> &parameters) override;
};

#endif // NEURALNETWORK_MONODENDRITESYMMETRICSTDP_H
