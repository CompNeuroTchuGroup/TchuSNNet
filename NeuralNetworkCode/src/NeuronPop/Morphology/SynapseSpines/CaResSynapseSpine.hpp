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
    cadouble calciumFree{};

    const double kinasesTotal{}; // If this is set to zero, there is no synapse
    const double phosphatasesTotal{};
    const double calmodulinTotal{};// We just multiply this by the logistic function to get the active ones
    const double neurograninTotal{};//THIS IS STILL NOT DONE

    double calmodulinActive{};
    double kinasesActive{};
    double kinasesInactive{};
    double phosphatasesActive{};
    double phosphatasesInactive{};

    double resourcesOldStep{};
    double resourcesAvailable{};

    CaResSynapseSpine(double kTot, double nTot, double camTot);
    ~CaResSynapseSpine() override = default;
    // End of step
    void RestoreSpine(); // This function MUST run right before diffusion
    // Profile methods
    std::vector<double> GetIndividualSynapticProfile() const override;
    std::string         GetIndividualSynapticProfileHeaderInfo() const override;
};

#endif