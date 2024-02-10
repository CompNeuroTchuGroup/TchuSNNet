//
// Created by Antoni Bertolin on 01.10.23
//
#ifndef _MASS_ACTION_CALCIUM_AND_RESOURCE_BASED_BRANCH_STRUCT_HPP
#define _MASS_ACTION_CALCIUM_AND_RESOURCE_BASED_BRANCH_STRUCT_HPP

#include "./Branch.hpp"
// #include "../BranchedMorphology.hpp"
#include "../../../GlobalFunctions.hpp"
#include "../SynapseSpines/MACRbpSynapseSpine.hpp"
#include <execution>

// include the spine class
struct MACRbPBranch : public Branch {

  // std::vector<std::vector<double>> waitingMatrix;//First axis is the position and time 0. Second axis is the time 1, 2, 3. Receives calcium of both
  // prespikes and postspikes
  // We swap vector positions (0,1;1,2;2,3), then set to false last vector. std::iter_swap or std::swap. Test speed in this.
  // Synapse access
  std::vector<MACRbPSynapseSpine> MACRbPspines; // VECTOR ACCOUNTS FOR EMPTY SYNAPSE SLOTS

  // Constants
  Constants constants;

  // int prespikeDelay; This is in Morhpo, given to constructor for matrix. No need to use afterwards.
  // Misc

  // Methods
  // Setup
  MACRbPBranch(std::vector<int> anteriorBranches, double gap, double branchLength, int branchId, Constants constants, MACRbPSynapseSpine spine);
  MACRbPBranch(double gap, double branchLength, int branchId, Constants constants, MACRbPSynapseSpine spine);
  void PostConnectSetUp(std::vector<BranchedSpinePtr> spineData) override;
  void PostSpikeCalciumFlux();
  // Input methods
  //  void PreSpikeCalciumInflux(TStepInt&& timestep);//We do swaps of the matrix and fill with zeroes. Matrix contains the calcium concentration
  //  already and we just add. void PreSpikeCalciumInflux();//We do swaps of the matrix and fill with zeroes. Matrix contains the calcium
  //  concentration already and we just add.
  // Reaction methods
  void Advect(); // All done in the scope of the branch!

  // Data retrieval
  double GetTotalWeight() const;
};
#endif