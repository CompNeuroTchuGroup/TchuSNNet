#include "CoopSynapseSpine.hpp"

// CoopSynapseSpine::CoopSynapseSpine(double distToSoma, double lastSpike, double weight):distToSoma{distToSoma}, BaseSynapseSpine(weight,lastSpike){
// }

std::valarray<double> CoopSynapseSpine::GetIndividualSynapticProfile() const
{
    std::valarray<double> dataArray(4);
    dataArray[0] = this->distToSoma;
    dataArray[1] = this->theta;
    dataArray[2] = this->weight;
    dataArray[3] = this->lastSpike;
    return dataArray;
}

std::string CoopSynapseSpine::GetIndividualSynapticProfileHeaderInfo() const
{
    return std::string("{<dist to soma>, <hetero cooperativity>, <weight>, <last spike>}");
}
