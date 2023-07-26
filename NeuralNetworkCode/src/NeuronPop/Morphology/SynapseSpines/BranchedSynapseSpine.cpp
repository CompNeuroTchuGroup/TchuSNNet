//
// Created by Antoni Bertolin on 14.06.23
//
#include "BranchedSynapseSpine.hpp"

// BranchedSynapseSpine::BranchedSynapseSpine(int distanceFromNode, double lastSpike, double weight, int branchId, int branchPositionId): distanceFromNode{ distanceFromNode }, branchId{ branchId }, branchPositionId {branchPositionId}, BaseSynapseSpine(weight, lastSpike)
// {
// }

std::vector<double> BranchedSynapseSpine::GetIndividualSynapticProfile() const
{
    std::vector<double> dataArray(3);
    dataArray.at(0) = this->distanceFromNode;
    dataArray.at(1) = this->branchId;
    dataArray.at(2) = this->weight;
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
