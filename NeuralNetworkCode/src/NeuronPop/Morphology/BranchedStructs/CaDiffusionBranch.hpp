//
// Created by Antoni Bertolin on 01.10.23
//
#ifndef _ALPHA_BRANCH_STRUCT_HEADER_
#define _ALPHA_BRANCH_STRUCT_HEADER_

#include "./Branch.hpp"
#include "../BranchedMorphology.hpp"
#include "../SynapseSpines/CaResSynapseSpine.hpp"
struct Constants{
    cadouble caDiffusionFct;//Consider delta x squared already here

    cadouble caDecayFct;//Consider already exponentiated
     
    double kinasesTotal;//Upper limit of Kinases (determines, with active, unactive species)
    double calcineurinTotal;//Upper limit of phosphatases (determines, with active, unactive species)
    double calmodulinTotal;// We just multiply this by the logistic function to get the active ones
    double neurograninTotal;//THIS IS STILL NOT DONE


    double reaction1Ctt;//k1, from K inactive to K active, using CaM
    double reaction2Ctt;//k2, from K active to K inactive, using CaM
    double reaction3Ctt;//k3, from N inactive to N active, using CaM
    double reaction4Ctt;//k4, from N active to N inactive, using CaM
    double reaction5Ctt;//kK, from resources to weight
    double reaction6Ctt;//kN, from weight to resources
    double reaction7Ctt;//From CaM and Ng to CaMNg
    double reaction8Ctt;//From CaMNg to CaM and Ng
    double reaction9Ctt;//From CaM and Ca to active CaM
    double reaction10Ctt;//From active CaM to CaM and Ca
    double reaction11Ctt;//From K bound to CaM to phosphorylated K
    double reaction12Ctt;//From phosphorylated K to K bound to CaM

    
    double resourceDiffusionFct;//Consider delta x squared already here
    // cadouble initialCalcium; 
    double initialResources, initialWeight;
};
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
    const double calcineurinTotal;//Upper limit of phosphatases (determines, with active, unactive species)
    const double calmodulinTotal;// We just multiply this by the logistic function to get the active ones
    const double neurograninTotal;//THIS IS STILL NOT DONE

    const double reaction1Ctt;//k1, from K inactive to K bound to CaM, using CaM
    const double reaction2Ctt;//k2, from K bound to CaM to K inactive, using CaM
    const double reaction3Ctt;//k3, from N inactive to N active, using CaM
    const double reaction4Ctt;//k4, from N active to N inactive, using CaM

    const double reaction5Ctt;//kK, from resources to weight
    const double reaction6Ctt;//kN, from weight to resources

    const double reaction7Ctt;//From CaM and Ng to CaMNg
    const double reaction8Ctt;//From CaMNg to CaM and Ng
    const double reaction9Ctt;//From CaM and Ca to active CaM
    const double reaction10Ctt;//From active CaM to CaM and Ca

    const double reaction11Ctt;//From K bound to CaM to phosphorylated K
    const double reaction12Ctt;//From phosphorylated K to K bound to CaM

    const double resourceDiffusionFct;//Consider delta x squared already here

    // int prespikeDelay; This is in Morhpo, given to constructor for matrix. No need to use afterwards.
    //Misc

    //Methods
    //Setup
    CaDiffusionBranch(std::vector<int>anteriorBranches,double gap, double branchLength, int branchId, TStepInt prespikeDelaySteps, Constants constants);
    CaDiffusionBranch(double gap, double branchLength, int branchId, TStepInt prespikeDelaySteps, Constants constants);
    void PostConnectSetUp(std::vector<BranchedSpinePtr> spineData) override;
    //Input methods
    void PreSpikeCalciumInflux(TStepInt&& timestep);//We do swaps of the matrix and fill with zeroes. Matrix contains the calcium concentration already and we just add.
    // void PreSpikeCalciumInflux();//We do swaps of the matrix and fill with zeroes. Matrix contains the calcium concentration already and we just add.    
    //Reaction methods
    void Advect(TStepInt step);//All done in the scope of the branch!
};
#endif