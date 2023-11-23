#include "./MACRbpSynapseSpine.hpp"

MACRbPSynapseSpine::MACRbPSynapseSpine(): connected{false} {
}

void MACRbPSynapseSpine::PreDiffusion() {
    //This has to run before reations, and before diffusion. Makes no sense to do both, and makes no sense to actually store the inactive forms in the spine, if we are going to calculate it to the timestep anyway
    // kinasesInactive=kinasesTotal-(kinasesCaM+kinasesPhospho);
    // phosphatasesInactive=calcineurinTotal-phosphatasesInactive;
    if(calciumFree<0 && resourcesAvailable<0){
        throw "Negative calcium and resources happened";
    }else if (calciumFree<0){
        throw "Negative calcium happened";
    } else if (resourcesAvailable<0){
        throw "Negative resources happened";
    }
    calciumOldStep=calciumFree;
    resourcesOldStep=resourcesAvailable;
}

std::vector<double> MACRbPSynapseSpine::GetIndividualSynapticProfile() const {
    std::vector<double> dataArray(5);
    dataArray.at(0) = this->branchId;
    dataArray.at(1) = this->weight;
    dataArray.at(2)= this->calciumFree;
    dataArray.at(3) = this->branchPositionId;
    dataArray.at(4) = this->prePopId;
    return dataArray;
}

std::string MACRbPSynapseSpine::GetIndividualSynapticProfileHeaderInfo() const {
    return std::string("{<branch ID>, <weight>, <calcium>, <position ID>, <presynaptic population ID>}");

}
