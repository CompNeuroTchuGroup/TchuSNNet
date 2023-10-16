#include "./ResourceCalciumDiffusionModel.hpp"#include "ResourceCalciumDiffusionModel.hpp"

BaseSpinePtr ResourceCalciumDiffusionModel::AllocateNewSynapse(const BranchTargeting &branchTarget) {
    int branch{AllocateBranch(branchTarget)};
    int position{PopSynapseSlotFromBranch(branch, branchTarget.firstSlotTrueLastSlotFalse)};
    // caDiffBranches.at(branch).CaDiffSpines.push_back(CaResSynapseSpine(kinasesTotal, calcineurinTotal, initialWeights));
    CaRsSpinePtr newSpine = &caDiffBranches.at(branch).CaResSpines.at(position);
    this->caResSpines.push_back(newSpine);
    // this->weightsSum += newSynapse->GetWeight();
    newSpine->idInMorpho=(this->baseSpineData.size());//this->spineIdGenerator++
    newSpine->weight=initialWeights;
    // Branch
    newSpine->branchPositionId=(position);
    newSpine->branchId=(branch);

    newSpine->connected=true;

    branches.at(branch)->synapseSlotClosedIndex.push_back(position);//Do we really need this?

    // Storage (other)
    this->baseSpineData.push_back(static_cast<BaseSpinePtr>(newSpine));
    this->branchedSpineData.push_back(static_cast<BranchedSpinePtr>(newSpine));

    return this->baseSpineData.back();
}
