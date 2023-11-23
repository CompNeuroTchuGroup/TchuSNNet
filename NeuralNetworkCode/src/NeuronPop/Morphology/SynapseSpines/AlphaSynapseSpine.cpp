//
// Created by Antoni Bertolin on 14.06.23
//
#include "AlphaSynapseSpine.hpp"



std::vector<double> AlphaSynapseSpine::GetIndividualSynapticProfile() const
{
    std::vector<double> dataArray(5);
    dataArray.at(0) = this->branchId;
    dataArray.at(1) = this->weight;
    dataArray.at(2) = this->alphaResources;
    dataArray.at(3) = this->branchPositionId;
    dataArray.at(4) = this->prePopId;
    return dataArray;
}

std::string AlphaSynapseSpine::GetIndividualSynapticProfileHeaderInfo() const{
    return std::string("{<branch ID>, <weight>, <alpha resources>, <position ID>, <presynaptic population ID>}");
}