//
// Created by Antoni Bertolin on 14.06.23
//
#ifndef _SYNAPSE_SPINE_BRANCHED_CLASS_HEADER_
#define _SYNAPSE_SPINE_BRANCHED_CLASS_HEADER_

#include "BaseSynapseSpine.hpp"

class BranchedSynapseSpine : public BaseSynapseSpine {

    protected:
    int branchId{}; //This has to be discrete
    int branchPositionId{}; //This has to be discrete
    double distanceFromNode{};//as the distance is id*gap, unnecessary
    //int uniqueTreeId{}; Useful only in the indexing of the triangular matrix if implemented
    bool isBranchedBool{true};

    public:
    BranchedSynapseSpine() = default;
    ~BranchedSynapseSpine() override = default;
    //BranchedSynapseSpine(int distanceFromNode, double lastSpike, double weight, int branchId, int branchPositionId);
    //getters
    int GetBranchId() const {return branchId;};
    int GetBranchPositionId() const {return branchPositionId;};
    double GetDistanceFromNode() const {return distanceFromNode;};
    bool GetBranchedBool() const { return isBranchedBool;}
    //setters
    void SetBranchId(int idIn){branchId=idIn;};
    void SetBranchPositionId(int positionId){branchPositionId=positionId;};
    void SetDistanceFromNode(double distance){distanceFromNode=distance;};

    //Profile function
    std::vector<double> GetIndividualSynapticProfile() const override;
    std::string GetIndividualSynapticProfileHeaderInfo() const override;
    //Friend functions
    //friend bool BranchIDCompare (const BranchedSpinePtr spine1, const BranchedSpinePtr spine2);
};

//bool BranchIDCompare (const BranchedSpinePtr spine1, const BranchedSpinePtr spine2);

#endif