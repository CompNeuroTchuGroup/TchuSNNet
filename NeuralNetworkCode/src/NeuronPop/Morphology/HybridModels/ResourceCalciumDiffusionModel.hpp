//
// Created by Antoni Bertolin on 14.06.23
//
#ifndef RESOURCE_CALCIUM_DIFFUSION_HETERO_SYNAPTIC_PLASTICITY_STDP_TRACE_BASED_HEADER_
#define RESOURCE_CALCIUM_DIFFUSION_HETERO_SYNAPTIC_PLASTICITY_STDP_TRACE_BASED_HEADER_

// List of forward declarations needed to break circular dependencies
//  class BranchedMorphology;
struct BranchTargeting;
#include "../../../GlobalFunctions.hpp"
#include "../BranchedMorphology.hpp"
#include "../BranchedStructs/CaDiffusionBranch.hpp"
#include "../SynapseSpines/CaResSynapseSpine.hpp"
#include <numeric>
using CaResSpinePtr = CaResSynapseSpine *;

class ResourceCalciumDiffusionModel : public BranchedMorphology {
    // This class models changes in weight depending on the concentration of kinases and phosphatases
  protected:
    std::vector<CaResSpinePtr>     caResSpines;
    std::vector<CaDiffusionBranch> caDiffBranches;
    // Constants
    double   kinasesTotal;         // Upper limit of Kinases (determines, with active, unactive species)
    double   phosphatasesTotal;    // Upper limit of phosphatases (determines, with active, unactive species)
    double   reaction1Ctt;         // k1, from K inactive to K active
    double   reaction2Ctt;         // k2, from K active to K inactive
    double   reaction3Ctt;         // k3, from N inactive to N active
    double   reaction4Ctt;         // k4, from N active to N inactive
    double   reaction5Ctt;         // kK, from resources to weight
    double   reaction6Ctt;         // kN, from weight to resources
    cadouble   caDiffusionFct;       // Consider delta x not squared yet here
    double   resourceDiffusionFct; // Consider delta x not squared yet here
    cadouble   caDecayFct;           // Consider not yet exponentiated
    cadouble initialCalcium;
    double   initialResources;
    cadouble prespikeCalcium, postspikeCalcium;
    TStepInt preSpikeDelaySteps;//either we calculate the modulus here or we keep to the swaps.
    TStepInt TStepMod;//Here we store the modulus (only calculated once). We will calculate it every timestep just for safety.

  public:
    // main Methods
    ResourceCalciumDiffusionModel() = default;
    explicit ResourceCalciumDiffusionModel(GlobalSimInfo *infoGlobal);
    ~ResourceCalciumDiffusionModel() override = default;

    void LoadParameters(const std::vector<FileEntry> &morphologyParameters) override;
    void CheckParameters(const std::vector<FileEntry> &parameters) override;
    void SaveParameters(std::ofstream &wParameterStream, std::string neuronIdentificator) const override;

    int  CreateBranch(std::vector<int> anteriorBranches) override;
    // void SetUpHashTable(); // Maybe if we want the delayed prespike to have a shape

    std::string GetType() const override { return IDstringTraceResourceHSTDP; };

    // Advect methods
    void Advect() override;
    // Pairing functions
    //  void UpdateCoopTrace(const ResourceTraceBranch* const branch);
    // bool CheckIfThereIsPairing(RBranchPtr branch, int synapseIDinBranch);
    void ApplyCoopTraceSpatialProfile(int branchSpineID, AlphaBranch &const branchID);
    // double CallKernelHashTable(int distanceToCenterInGaps);
    // Plasticity events functions
    // void ApplyEffects();//Here we increase the plasticity count of synapse and branch
    // void STDPPotentiation(ResourceSpinePtr synapse);
    // void STDPDepression(ResourceSpinePtr synapse);
    // Reset methods
    void Reset() override; // Wrapper plus clearing some of the vectors. Last Reset method to run in chronological
                           // order, where we call the ticks and the general upkeep
    // void DeleteEffects();//Here, if counter==countMax, erase in that index the element of every vector (first store
    // index, then REVERSE remove the removelist indexes with .rbegin and .rend instead of .begin and .end)
    // Container should be ordered by definition, but std::sort(array.begin(), array.end()) would ensure so.
    void DecayAllTraces(); // Last method called in Reset()
    // void ClearSynapseSets();
    // Recalc methods. These methods have to be done per branch
    void ComputeAlphas(AlphaBranch &const branch);    // Run in LP
    void ComputeWeights(AlphaBranch &const branch);   // Run in LP
    void ComputeAlphaSums(AlphaBranch &const branch); // Called inside recalc weights
    // Record methods
    void RecordPostSpike() override;
    void RecordExcitatoryPreSpike(int spikedSpineId) override; // Here set the trigger count to 0
    void PostConnectSetUp() override;
    // Allocation methods
    BaseSpinePtr AllocateNewSynapse(
        const BranchTargeting &branchTarget) override; // Call the Branched one inside before setting all counters
                                                  // Remember to set all counts to maxCount
    // Record functions
    std::vector<double> GetOverallSynapticProfile() const override;
    std::string         GetOverallSynapticProfileHeaderInfo() const override;
    // void CalcMorphoPlasticityEvents() override;
    // For debugging purposes
    bool IgnoreJDParameters() const override { return true; }
};

#endif