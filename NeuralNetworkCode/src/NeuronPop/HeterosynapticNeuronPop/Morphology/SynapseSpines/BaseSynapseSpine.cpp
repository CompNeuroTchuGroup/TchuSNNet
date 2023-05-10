#include "BaseSynapseSpine.hpp"


// BaseSynapseSpine::BaseSynapseSpine(double weight, double lastSpike): weight{weight}, lastSpike{lastSpike}
// {
// }

bool MorphoIDCompare(const BaseSpinePtr &spine1, const BaseSpinePtr &spine2)
{
    return (spine1->idInMorpho<spine2->idInMorpho);
}
