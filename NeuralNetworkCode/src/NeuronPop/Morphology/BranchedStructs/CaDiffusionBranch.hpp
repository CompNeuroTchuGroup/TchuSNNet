//
// Created by Antoni Bertolin on 01.10.23
//
#ifndef _ALPHA_BRANCH_STRUCT_HEADER_
#define _ALPHA_BRANCH_STRUCT_HEADER_

#include "./Branch.hpp"
#include "../BranchedMorphology.hpp"
#include "../SynapseSpines/CaResSynapseSpine.hpp"
// include the spine class
struct CaDiffusionBranch : public Branch {

    std::vector<std::vector<cadouble>> waitingMatrix;//First axis is the position and time 0. Second axis is the time 1, 2, 3. Receives calcium of both prespikes and postspikes
    //We swap vector positions (0,1;1,2;2,3), then set to false last vector. std::iter_swap or std::swap. Test speed in this.
    //Synapse access
    std::vector<CaResSynapseSpine> CaResSpines;//VECTOR ACCOUNTS FOR EMPTY SYNAPSE SLOTS

    //Constants
    const cadouble caDiffusionFct;//Consider delta x squared already here

    const cadouble caDecayFct;//Consider already exponentiated
     
    const double kinasesTotal;//Upper limit of Kinases (determines, with active, unactive species)
    const double phosphatasesTotal;//Upper limit of phosphatases (determines, with active, unactive species)
    const double reaction1Ctt;//k1, from K inactive to K active, using CaM
    const double reaction2Ctt;//k2, from K active to K inactive, using CaM
    const double reaction3Ctt;//k3, from N inactive to N active, using CaM
    const double reaction4Ctt;//k4, from N active to N inactive, using CaM
    const double camSaturationSteep;//This constant is for the CaM saturation curve.
    const double camSaturationInflection;//This constant is for the CaM saturation curve
    const double reaction5Ctt;//kK, from resources to weight
    const double reaction6Ctt;//kN, from weight to resources
    const double resourceDiffusionFct;//Consider delta x squared already here

    // int prespikeDelay; This is in Morhpo, given to constructor for matrix. No need to use afterwards.
    //Misc

    //Methods
    //Setup
    CaDiffusionBranch(std::vector<int>anteriorBranches,double gap, double branchLength, int branchId, cadouble initialCalcium, double initialResources, double initialWeight, double kinasesTotal, double phosphatasesTotal, double reaction1Ctt, double reaction2Ctt, double reaction3Ctt, double reaction4Ctt, double reaction5Ctt, double reaction6Ctt,cadouble caDiffusionFct,double resourceDiffusionFct,cadouble caDecayFct, int prespikeDelaySteps,double camSaturationSteep, double camSaturationInflection);
    CaDiffusionBranch(double gap, double branchLength, int branchId, cadouble initialCalcium, double initialResources, double initialWeight, double kinasesTotal, double phosphatasesTotal, double reaction1Ctt, double reaction2Ctt, double reaction3Ctt, double reaction4Ctt, double reaction5Ctt, double reaction6Ctt,cadouble caDiffusionFct,double resourceDiffusionFct,cadouble caDecayFct, int prespikeDelaySteps,double camSaturationSteep, double camSaturationInflection);    
    void PostConnectSetUp(std::vector<BranchedSpinePtr> spineData) override;
    //Input methods
    void PreSpikeCalciumInflux(TStepInt timestep);//We do swaps of the matrix and fill with zeroes. Matrix contains the calcium concentration already and we just add.
    // void PreSpikeCalciumInflux();//We do swaps of the matrix and fill with zeroes. Matrix contains the calcium concentration already and we just add.    
    //Reaction methods
    void ExecuteReactions();//All done in the scope of the branch!
    void CalciumCheck();
};
#endif