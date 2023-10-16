//
// Created by Antoni Bertolin on 14.06.23
//
#ifndef _RESOURCE_SYNAPSE_SPINE_CLASS_HEADER_
#define _RESOURCE_SYNAPSE_SPINE_CLASS_HEADER_

#include "BranchedSynapseSpine.hpp"
using cadouble = double; //If we go to long for precision concerns

struct CaResSynapseSpine : public BranchedSynapseSpine {

  public:
    bool connected{false};

    cadouble calciumOldStep{};
    cadouble calciumFree{0.08};

    double calmodulinActive{};
    double calmodulinNeurogranin{};

    double kinasesCaM{};
    double kinasesPhospho{};
    // double kinasesInactive{};

    double calcineurinActive{};
    // double phosphatasesInactive{};

    double resourcesOldStep{};
    double resourcesAvailable{};

    CaResSynapseSpine();
    ~CaResSynapseSpine() override = default;
    // End of step
    void PreDiffusion(); // This function MUST run right before diffusion
    // Profile methods
    std::vector<double> GetIndividualSynapticProfile() const override;
    std::string         GetIndividualSynapticProfileHeaderInfo() const override;
};

#endif