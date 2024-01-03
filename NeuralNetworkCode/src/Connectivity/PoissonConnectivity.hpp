#ifndef _POISSON_CONNECTIVITY_HPP
#define _POISSON_CONNECTIVITY_HPP
#include "../GlobalFunctions.hpp"
#include "Connectivity.hpp"
#include <iostream>
#include <random>
#include <vector>

class Synapse;

class PoissonConnectivity : public Connectivity {
protected:
  signed long totalConnections{};
  double      connectionLambda{};

public:
  PoissonConnectivity(Synapse *synapse, GlobalSimInfo *infoGlobal);
  ~PoissonConnectivity() override = default;

  void        ConnectNeurons() override;
  std::string GetTypeStr() const override { return IDstringPoissonConnectivity; }

  double GetExpectedConnections() const override;

  void SaveParameters(std::ofstream &wParameterStream, std::string idString) const override;
  void LoadParameters(const std::vector<FileEntry> &parameters) override;
};

#endif /* PoissonConnectivity_hpp */
