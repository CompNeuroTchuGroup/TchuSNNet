//
// Created by Antoni Bertolin on 14.06.23
//
#ifndef MASS_ACTION_CALCIUM_AND_RESOURCE_BASED_SYNAPSE_SPINE_CLASS_HEADER_
#define MASS_ACTION_CALCIUM_AND_RESOURCE_BASED_SYNAPSE_SPINE_CLASS_HEADER_

#include "BranchedSynapseSpine.hpp"
// using double = double; //If we go to long for precision concerns

struct MACRbPSynapseSpine : public BranchedSynapseSpine {

  public:
    bool connected{false};

    double preTransientIncrease{};//Abstract trace
    double preTransient{};//Calcium trace (interacts with postspike)

    double postTransientIncrease{};//Abstract trace
    double postTransient{};//Calcium trace (interacts with postspike)

    double calciumOldStep{};
    double calciumFree{};

    double calmodulinActive{};
    double calmodulinNeurogranin{};

    double kinasesCaM{};
    double kinasesPhospho{};
    // double kinasesInactive{};

    double calcineurinActive{};
    // double phosphatasesInactive{};

    double resourcesOldStep{};
    double resourcesAvailable{};

    MACRbPSynapseSpine();
    ~MACRbPSynapseSpine() override = default;
    // End of step
    void PreDiffusion(); // This function MUST run right before diffusion
    // Profile methods
    std::vector<double> GetIndividualSynapticProfile() const override;
    std::string         GetIndividualSynapticProfileHeaderInfo() const override;
};

#endif