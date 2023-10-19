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
    std::vector<double> dataArray(5);
    dataArray.at(0) = this->branchId;
    dataArray.at(1) = this->weight;
    dataArray.at(2)= this->calciumFree;
    dataArray.at(3) = this->branchPositionId;
    dataArray.at(4) = this->prePopId;
    return dataArray;
}

std::string CaResSynapseSpine::GetIndividualSynapticProfileHeaderInfo() const {
    return std::string("{<branch ID>, <weight>, <calcium>, <position ID>, <presynaptic population ID>}");

}
