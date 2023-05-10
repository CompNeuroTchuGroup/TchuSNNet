#include "BranchedSynapseSpine.hpp"

// BranchedSynapseSpine::BranchedSynapseSpine(int distanceFromNode, double lastSpike, double weight, int branchId, int branchPositionId): distanceFromNode{ distanceFromNode }, branchId{ branchId }, branchPositionId {branchPositionId}, BaseSynapseSpine(weight, lastSpike)
// {
// }

std::valarray<double> BranchedSynapseSpine::GetIndividualSynapticProfile() const
{
    std::valarray<double> dataArray(3);
    dataArray[0] = this->distanceFromNode;
    dataArray[1] = this->branchId;
    dataArray[2] = this->weight;
//    dataArray[3] = this->lastSpike;
    return dataArray;
}

std::string BranchedSynapseSpine::GetIndividualSynapticProfileHeaderInfo() const
{
    return std::string("{<dist to branch node>, <branch ID>, <weight>, <last spike>}");
}

// bool BranchIDCompare(const BranchedSpinePtr &spine1, const BranchedSpinePtr &spine2)
// {
//     return (spine1->branchPositionId<spine2->branchPositionId);
// }
