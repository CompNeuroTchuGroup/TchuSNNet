//
// Created by Antoni Bertolin on 14.06.23
//
#ifndef BRANCHED_MORPHOLOGY_HEADER_
#define BRANCHED_MORPHOLOGY_HEADER_
//List of forward declarations needed to break circular dependencies
// class Morphology;
// class Branch;
// class BranchTargeting;
#include "./Morphology.hpp"
#include "./SynapseSpines/BaseSynapseSpine.hpp"
#include "./SynapseSpines/BranchedSynapseSpine.hpp"
#include "./BranchedStructs.hpp"
#include <string>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <random>


class BranchedMorphology : public Morphology {

protected:
//MonoDendriteSTDP moved
    signed long spineIdGenerator{0}; // variable used to allocate new synapses. type is legacy compatible
    int branchIdGenerator{0}; //Same but in branches. Should be 1 or 0?
    //std::vector<bool> integratePostSpike;// Not necessary
    //std::vector<bool> integratePreSpike;//Not necessary

    //Weight normalization vars
    double synapticGap{};
    double branchLength{};

//Branched specific
    int branchings{1};

    bool orderedSpineAllocationB{false};// If not properly loaded from LP, exception
    bool randomSpineAllocationB{false};


    /*bool setBranchAllocationB{false};
    bool OrderedBranchAllocationB{false};// If not properly loaded from LP, exception
    bool RandomBranchAllocationB{false};*/

    std::vector<BranchPtr> branches{};//unique_ptr's constructor is explicit, so you either need to use emplace_back or stuff.push_back(std::unique_ptr<int>(new static_cast<int>(i)));. Between the two, emplace_back is much cleaner.
    std::vector<BranchedSpinePtr> branchedSpineData;//They are just pointers, what is the worst that can happen by having multiple copies?
public:
    explicit BranchedMorphology(GlobalSimInfo* infoGlobal);
    ~BranchedMorphology() = default;

    //Methods derived from the MonoDendriteSTDP and Morphology classes
    void SaveParameters(std::ofstream& wParameterFile, std::string neuronIdentificator) const override;
    void LoadParameters(const std::vector<FileEntry>& parameters) override;
    void CheckParameters(const std::vector<FileEntry>& parameters) override;
    virtual std::string GetType() const = 0;

    void RecordPostSpike() override;
    void RecordExcitatoryPreSpike(int spikedSynapseId) override; 
    //virtual std::vector<double> GetOverallSynapticProfile() const;
    //void CalcMorphoPlasticityEvents() override;
    
    virtual void Advect()=0;
    void Reset() override;

    //Branched specific methods
    void SetUpBranchedMorphology();
    void SetUpSynapseSlots(BranchPtr branch); //This function will set up the open synapse slots of a branch object with its id.This one I have to define in the parallel synaptic connectivity masks or the derived classes
    //setUp SYnapse slots is called for every branch in a loop and depending on the bool (universal for all branches for now) it calls random or ordered.
    //The overriding function calls functions of BMorpho. 
    virtual void SetUpBranchings(int remainingBranchingEvents, std::vector<int> anteriorBranches = std::vector<int>());// Here we set up the vector with the branches
    void PostConnectSetUp() override;
    //Allocation shennanigans
    int GetNoBranches(){return static_cast<int>(branches.size());}
    BaseSpinePtr AllocateNewSynapse(const BranchTargeting& synapse) override = 0; //Use the reference to call GetBranchTarget

    int AllocateBranch(const BranchTargeting synapse);//The selected branch allocation is simple. This function is called in AllocateNewSynapse
    int RandomBranchAllocation();
    int OrderedBranchAllocation();
    //setBranchAllocation() is implicit in the function (or has to be) allocate NewSynapse
    //virtual int orderedGuidedBranchAllocation(const char DendriticSubRegionID);
    void RandomSynapseAllocation(BranchPtr branch);
    void OrderedSynapseAllocation(BranchPtr branch);//These two are coming from the SetUpSynapseSlots already, called depending on a bool. 
    int PopSynapseSlotFromBranch(int branch, bool firstSlotTrueLastSlotFalse);
    //virtual void AlternatedSynapseAllocation(BranchPtr branch);
    //
    int GenerateBranchId(){return branchIdGenerator++;}

};

#endif