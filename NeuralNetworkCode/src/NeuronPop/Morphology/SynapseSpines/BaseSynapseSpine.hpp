//
// Created by Antoni Bertolin on 14.06.23
//
#ifndef _BASE_SYNAPSE_SPINE_CLASS_HPP
#define _BASE_SYNAPSE_SPINE_CLASS_HPP

#include "../../../GlobalFunctions.hpp"
#include <string>
#include <vector>

struct BaseSynapseSpine {

  // Legacy variables
  signed long idInMorpho{}; // id for synapse within its population

  PopInt prePopId{};
  double couplingStrength{1};
  double weight{}; // The negative weight comes from J, weight is just a factor to multiply (for now)
  // double lastSpike{};

  // Branched variables
  // Synapse related

  // Necessary mutex

  // Constructors
  BaseSynapseSpine() = default;
  BaseSynapseSpine(double weight);
  virtual ~BaseSynapseSpine() = default;
  // BaseSynapseSpine(double weight, double lastSpike);// Not currently in use
  // Methods
  // Getters
  double GetWeight() const { return weight * couplingStrength; };
  double GetWeightUncoupled() const { return weight; };
  // Recorder functions
  virtual std::vector<double> GetIndividualSynapticProfile() const           = 0;
  virtual std::string         GetIndividualSynapticProfileHeaderInfo() const = 0;
};

#endif