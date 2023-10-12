#include "CaResSynapseSpine.hpp"

CaResSynapseSpine::CaResSynapseSpine(double kTot, double nTot, double camTot): kinasesTotal{kTot}, phosphatasesTotal{nTot},calmodulinTotal{camTot}, connected{false} {
}

void CaResSynapseSpine::RestoreSpine() {
    kinasesInactive=kinasesTotal-kinasesActive;
    phosphatasesInactive=phosphatasesTotal-phosphatasesInactive;
    calciumOldStep=calciumFree;
    resourcesOldStep=resourcesAvailable;
}

std::vector<double> CaResSynapseSpine::GetIndividualSynapticProfile() const {
    return std::vector<double>();
}

std::string CaResSynapseSpine::GetIndividualSynapticProfileHeaderInfo() const {
    return std::string();
}
