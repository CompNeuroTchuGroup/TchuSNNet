//
// Created by Antoni Bertolin on 14.06.23
//
#include "CoopSynapseSpine.hpp"

// CoopSynapseSpine::CoopSynapseSpine(double distToSoma, double lastSpike, double weight):distToSoma{distToSoma}, BaseSynapseSpine(weight,lastSpike){
// }

std::vector<double> CoopSynapseSpine::GetIndividualSynapticProfile() const {
  std::vector<double> dataArray(4);
  dataArray.at(0) = this->distToSoma;
  dataArray.at(1) = this->theta;
  dataArray.at(2) = this->weight;
  dataArray.at(3) = this->lastSpike;
  return dataArray;
}

std::string CoopSynapseSpine::GetIndividualSynapticProfileHeaderInfo() const {
  return std::string("{<dist to soma>, <hetero cooperativity>, <weight>, <last spike>}");
}
