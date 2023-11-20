//
// Created by Antoni Bertolin on 01.10.23
//
#ifndef _CALCIUM_RESOURCE_DIFFUSION_BRANCH_STRUCT_HEADER_
#define _CALCIUM_RESOURCE_DIFFUSION_BRANCH_STRUCT_HEADER_

#include "./Branch.hpp"
// #include "../BranchedMorphology.hpp"
#include "../../../GlobalFunctions.hpp"
#include "../SynapseSpines/CaResSynapseSpine.hpp"
struct Constants{ //currently 19 + the prespike/postspike calcium and the prespike calcium influx delay, so 21
    double caDiffusionFct;//Consider delta x squared already here

    double calciumExtrusionCtt;//Consider already exponentiated
     
    double kinasesTotal;//Upper limit of Kinases (determines, with active, unactive species)
    double calcineurinTotal;//Upper limit of phosphatases (determines, with active, unactive species)
    double calmodulinTotal;// We just multiply this by the logistic function to get the active ones
    double neurograninTotal;//THIS IS STILL NOT DONE

    double reaction1Ctt;//From CaM and Ng to CaMNg
    double reaction2Ctt;//From CaMNg to CaM and Ng
    double reaction3Ctt;//From CaM and Ca to active CaM
    double reaction4Ctt;//From active CaM to CaM and Ca
    double reaction5Ctt;//From N inactive to N active, using CaM
    double reaction6Ctt;//From N active to N inactive, using CaM
    double reaction7Ctt;//From K bound to CaM to phosphorylated K
    double reaction8Ctt;//From phosphorylated K to K bound to CaM
    double reaction9Ctt;//From K inactive to K active, using CaM
    double reaction10Ctt;//From K active to K inactive, using CaM
    double reaction11Ctt;//From resources to weight
    double reaction12Ctt;//From weight to resources

    
    double resourceDiffusionFct;//Consider delta x squared already here
    double resourceConversionFct;
    double calciumInfluxBasal;
    double initialResources, initialWeight;

    double preCalciumFluxFactor;
    double preCalciumRiseRate;
    double preCalciumDecayRate;

    double postCalciumFluxFactor;
    double postCalciumRiseRate;
    double postCalciumDecayRate;

    double nonlinearFactorNMDA;
};
// include the spine class
struct CaDiffusionBranch : public Branch {

    // std::vector<std::vector<double>> waitingMatrix;//First axis is the position and time 0. Second axis is the time 1, 2, 3. Receives calcium of both prespikes and postspikes
    //We swap vector positions (0,1;1,2;2,3), then set to false last vector. std::iter_swap or std::swap. Test speed in this.
    //Synapse access
    std::vector<CaResSynapseSpine> CaResSpines;//VECTOR ACCOUNTS FOR EMPTY SYNAPSE SLOTS

    //Constants
    Constants constants;

    // int prespikeDelay; This is in Morhpo, given to constructor for matrix. No need to use afterwards.
    //Misc

    //Methods
    //Setup
    CaDiffusionBranch(std::vector<int>anteriorBranches,double gap, double branchLength, int branchId, Constants constants);
    CaDiffusionBranch(double gap, double branchLength, int branchId, Constants constants);
    void PostConnectSetUp(std::vector<BranchedSpinePtr> spineData) override;
    void PostSpikeCalciumFlux();
    //Input methods
    // void PreSpikeCalciumInflux(TStepInt&& timestep);//We do swaps of the matrix and fill with zeroes. Matrix contains the calcium concentration already and we just add.
    // void PreSpikeCalciumInflux();//We do swaps of the matrix and fill with zeroes. Matrix contains the calcium concentration already and we just add.    
    //Reaction methods
    void Advect();//All done in the scope of the branch!
};
#endif