//
// Created by Antoni Bertolin on 14.06.23
//
#ifndef _RESOURCE_SYNAPSE_SPINE_CLASS_HPP
#define _RESOURCE_SYNAPSE_SPINE_CLASS_HPP

#include "BranchedSynapseSpine.hpp"
#include <string>
#include <vector>

struct AlphaSynapseSpine : public BranchedSynapseSpine {

  // Central vars
  double alphaResources{};

  // Stimulus vector methods
  //  void TickStimulusCounts();//Called in Reset(), but both should be mutually exclusive with AATE above
  //  void CullStimulusVectors();//Called in Reset()
  // Profile methods
  std::vector<double> GetIndividualSynapticProfile() const override;
  std::string         GetIndividualSynapticProfileHeaderInfo() const override;
};

#endif