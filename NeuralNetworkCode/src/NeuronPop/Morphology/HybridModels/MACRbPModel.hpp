//
// Created by Antoni Bertolin on 14.06.23
//
#ifndef _MASS_ACTION_CALCIUM_AND_RESOURCE_BASED_PLASTICITY_HPP
#define _MASS_ACTION_CALCIUM_AND_RESOURCE_BASED_PLASTICITY_HPP

// List of forward declarations needed to break circular dependencies
//  class BranchedMorphology;
struct BranchTargeting;

#include "../../../GlobalFunctions.hpp"
#include "../BranchedMorphology.hpp"
#include "../BranchedStructs/MACRbPBranch.hpp"
#include "../SynapseSpines/MACRbPSynapseSpine.hpp"
#include <numeric>

class MACRbPModel : public BranchedMorphology {
  // This will be called the Mass Action Calcium and Resource based Plasticity Model, MACRbP model
protected:
  Constants constants{};

  std::vector<MACRbpSpinePtr> caResSpines;
  std::vector<MACRbPBranch>   caDiffBranches;
  // Constants
  double prespikeCalcium{}, postspikeCalcium{};

  double preCalciumRiseTau{};
  double preCalciumDecayTau{};

  double postCalciumRiseTau{};
  double postCalciumDecayTau{};

  double availResourcesRatio{};

  double calciumBasal{};
  // TStepInt preSpikeDelaySteps{}; // either we calculate the modulus here
  // TStepInt TStepModded{};           // Here we store the modulus (only calculated once). We will calculate it every timestep just for safety.
  // TStepInt TStepInput{};//Here is where the input happens
public:
  // main Methods
  MACRbPModel() = default;
  explicit MACRbPModel(GlobalSimInfo *infoGlobal);
  ~MACRbPModel() override = default;

  void LoadParameters(const std::vector<FileEntry> &morphologyParameters) override;
  void CheckParameters(const std::vector<FileEntry> &parameters) override;
  void SaveParameters(std::ofstream &wParameterStream, std::string neuronIdentificator) const override;

  int CreateBranch(std::vector<int> anteriorBranches) override;

  std::string GetType() const override { return IDstringTraceResourceHSTDP; };

  // Advect methods
  void Advect() override;
  // Record methods
  void RecordPostSpike() override;
  void RecordExcitatoryPreSpike(BaseSpinePtr spinePtr) override; // Here set the trigger count to 0
  void PostConnectSetUp() override;
  // Allocation methods
  BaseSpinePtr AllocateNewSynapse(BranchTargeting &branchTarget) override; // Call the Branched one inside before setting all counters
  // Remember to set all counts to maxCount
  // Record functions
  std::vector<double> GetOverallSynapticProfile() const override;
  std::string         GetOverallSynapticProfileHeaderInfo() const override;
  // void CalcMorphoPlasticityEvents() override;
  // For debugging purposes
  bool IgnoreJDParameters() const override { return true; }
};

#endif