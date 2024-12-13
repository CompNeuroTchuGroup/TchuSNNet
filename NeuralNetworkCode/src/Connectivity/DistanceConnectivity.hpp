#ifndef _DISTANCE_CONNECTIVITY_HPP
#define _DISTANCE_CONNECTIVITY_HPP

#include "../GlobalFunctions.hpp"
#include "../Synapse/Synapse.hpp"
#include "Connectivity.hpp"
#include <iostream>
#include <random>
#include <vector>

class Synapse;

class DistanceConnectivity : public Connectivity {
protected:
  double peakProbability{};
  double stdProbability{1};
  double projectionLengthFactor{1}; //==a
  int    isConnectionExact{};

public:
  DistanceConnectivity(Synapse *synapse, GlobalSimInfo *infoGlobal);
  ~DistanceConnectivity() override = default;

  void        ConnectNeurons();
  void        ConnectNeuronsExact();
  double      GetExpectedConnections() const override;
  std::string GetTypeStr() const override { return IDstringDistanceConnectivity; }

  void SaveParameters(std::ofstream &wParameterStream, std::string idString) const;
  void LoadParameters(const std::vector<FileEntry> &parameters);
};

#endif /* DistanceConnectivity_hpp */
