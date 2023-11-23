//
// Created by Antoni Bertolin on 14.06.23
//
#ifndef _RESOURCE_SYNAPSE_SPINE_CLASS_HEADER_
#define _RESOURCE_SYNAPSE_SPINE_CLASS_HEADER_

#include "BranchedSynapseSpine.hpp"
#include <vector>
#include <list>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <unordered_map>

struct AlphaSynapseSpine : public BranchedSynapseSpine {

    //Central vars
    double alphaResources{};

    //Stimulus vector methods
    // void TickStimulusCounts();//Called in Reset(), but both should be mutually exclusive with AATE above
    // void CullStimulusVectors();//Called in Reset()
    //Profile methods
    std::vector<double> GetIndividualSynapticProfile() const override;
    std::string GetIndividualSynapticProfileHeaderInfo() const override;
};

#endif