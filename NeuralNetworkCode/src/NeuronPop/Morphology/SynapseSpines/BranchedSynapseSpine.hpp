//
// Created by Antoni Bertolin on 14.06.23
//
#ifndef _SYNAPSE_SPINE_BRANCHED_CLASS_HPP
#define _SYNAPSE_SPINE_BRANCHED_CLASS_HPP

#include "BaseSynapseSpine.hpp"
#include <string>
#include <vector>

struct BranchedSynapseSpine : public BaseSynapseSpine {

  int branchId{};         // This has to be discrete
  int branchPositionId{}; // This has to be discrete
  // int uniqueTreeId{}; Useful only in the indexing of the triangular matrix if implemented

  BranchedSynapseSpine() = default;
  // BranchedSynapseSpine(double weight);
  ~BranchedSynapseSpine() override = default;
  // BranchedSynapseSpine(int distanceFromNode, double lastSpike, double weight, int branchId, int branchPositionId);

  // Profile function
  std::vector<double> GetIndividualSynapticProfile() const override           = 0;
  std::string         GetIndividualSynapticProfileHeaderInfo() const override = 0;
  // Friend functions
  // friend bool BranchIDCompare (const BranchedSpinePtr spine1, const BranchedSpinePtr spine2);
};

// bool BranchIDCompare (const BranchedSpinePtr spine1, const BranchedSpinePtr spine2);

#endif