//
// Created by Antoni Bertolin on 14.06.23
//
#ifndef _SYNAPSE_SPINE_COOP_CLASS_HPP
#define _SYNAPSE_SPINE_COOP_CLASS_HPP

#include "BaseSynapseSpine.hpp"
#include <string>
#include <vector>

class CoopSynapseSpine : public BaseSynapseSpine {

protected:
  double distToSoma{};
  double theta{}; // heterosynaptic cooperativity
  double lastSpike{};

public:
  CoopSynapseSpine() = default;
  // CoopSynapseSpine(double distToSoma, double lastSpike, double weight);

  // getters
  double GetDistToSoma() const { return distToSoma; };
  double GetTheta() const { return theta; };
  double GetLastSpike() const { return lastSpike; };
  // setters
  void SetDistToSoma(double dist) { distToSoma = dist; };
  void SetTheta(double thetaIn) { theta = thetaIn; };
  void SetLastSpike(double time) { lastSpike = time; };
  // Misc
  void AddToTheta(double hEffect) { theta += hEffect; }

  // Profile function
  std::vector<double> GetIndividualSynapticProfile() const override;
  std::string         GetIndividualSynapticProfileHeaderInfo() const override;
};

#endif