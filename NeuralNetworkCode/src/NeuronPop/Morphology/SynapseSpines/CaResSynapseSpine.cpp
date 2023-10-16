#include "CaResSynapseSpine.hpp"

CaResSynapseSpine::CaResSynapseSpine(): connected{false} {
}

void CaResSynapseSpine::PreDiffusion() {
    //This has to run before reations, and before diffusion. Makes no sense to do both, and makes no sense to actually store the inactive forms in the spine, if we are going to calculate it to the timestep anyway
    // kinasesInactive=kinasesTotal-(kinasesCaM+kinasesPhospho);
    // phosphatasesInactive=calcineurinTotal-phosphatasesInactive;
    if (calciumFree<0 || resourcesAvailable<0){
        throw "Negative something happened";
    }
    calciumOldStep=calciumFree;
    resourcesOldStep=resourcesAvailable;
}

std::vector<double> CaResSynapseSpine::GetIndividualSynapticProfile() const {
    return std::vector<double>();
}

std::string CaResSynapseSpine::GetIndividualSynapticProfileHeaderInfo() const {
    return std::string();
}
