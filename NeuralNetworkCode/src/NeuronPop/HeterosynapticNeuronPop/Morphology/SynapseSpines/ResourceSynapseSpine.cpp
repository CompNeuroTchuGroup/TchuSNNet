#include "ResourceSynapseSpine.hpp"

std::valarray<double> ResourceSynapseSpine::GetIndividualSynapticProfile() const
{
    std::valarray<double> dataArray(4);
    dataArray[0] = this->branchId;
    dataArray[1] = this->weight;
    dataArray[2] = this->alphaResources;
    dataArray[3] = this->distanceFromNode;
    return dataArray;
}

std::string ResourceSynapseSpine::GetIndividualSynapticProfileHeaderInfo() const
{
    return std::string("{<branch ID>, <weight>, <alpha resources>, <distance to node>}");
}

void ResourceSynapseSpine::RecalculateAlphaResources()
{
    //alphaResources=(kBasal+kStimmulus)/(nBasal+nStimmulus);
    alphaResources=alphaBasal+alphaStimmulus;
    DecayAlphaResources();
}

void ResourceSynapseSpine::DecayAlphaResources()
{
    // kStimmulus*=kStimmulusExpDecay;
    // nStimmulus*=nStimmulusExpDecay;
    alphaStimmulus *= alphaStimmulusExpDecay;
}

void ResourceSynapseSpine::RecalcWeight(double weightResourceFactor)
{
    this->weight=this->alphaResources*weightResourceFactor;
}

void ResourceSynapseSpine::AddTempResourcesToSpine(double alphaStimmulusInput)
{
    // kStimmulusTempEffect.push_back(kStimmulusInput);
    // nStimmulusTempEffect.push_back(nStimmulusInput);
    //stimmulusEffectCount.push_back(0);
    depressionAlphaTempAndCount.push_back(alphaStimmulusInput);
    potentiationAlphaTempAndCount.emplace_back(PairDI(std::abs(alphaStimmulusInput), 0));
}

bool ResourceSynapseSpine::ApplyAllTemporaryEffectsOnPostspike(const DHashMap& STDPdecayMap)
{
    // kStimmulus=std::accumulate(kStimmulusTempEffect.begin(), kStimmulusTempEffect.end(), kStimmulus, [PotentiationDepressionRatio](double accumulator, double kStemp){return accumulator + kStemp*PotentiationDepressionRatio;});
    // nStimmulus=std::accumulate(nStimmulusTempEffect.begin(), nStimmulusTempEffect.end(), nStimmulus, [PotentiationDepressionRatio](double accumulator, double nStemp){return accumulator + nStemp*PotentiationDepressionRatio;});
    // kStimmulusTempEffect.clear();
    // nStimmulusTempEffect.clear();
    //stimmulusEffectCount.clear();
    //PotentiationDepressionRatio is to boost potentiation compared to depression, but also to input the current decay
    //The next line goes over all the pairs in the potentiation vector and adds their stored alpha stimmulus (stored there by pairing events) and adds them to the alphaStimmulus variable
    //In addition, it multiplies with the decay to postspike of those effects and multiplies by the potentiation/depression ratio (a parameter from .txt)
    alphaStimmulus=std::accumulate(potentiationAlphaTempAndCount.begin(), potentiationAlphaTempAndCount.end(), alphaStimmulus, [this, STDPdecayMap](double accumulator, const PairDI& alphaStemp){return accumulator + alphaStemp.first*PotentiationDepressionRatio*STDPdecayMap.at(alphaStemp.second);});
    // if (alphaStimmulus+alphaBasal<0.0){
    //     alphaStimmulus= (-alphaBasal);
    // }
    //The bool return is to keep track of plasticity events (for now not very useful, as still no distinction by branch)
    return (!potentiationAlphaTempAndCount.size()==0);
}

// void ResourceSynapseSpine::ApplyAllTempEffectsOnConflictPotentiation(double PotentiationDepressionRatio)
// {
//     //CUrrently I have no way of properly calculating the decay of potentiation proper by calculating the area under the graph.
//     alphaStimmulus=std::accumulate(potentiationAlphaTempAndCount.begin(), potentiationAlphaTempAndCount.end(), alphaStimmulus, [PotentiationDepressionRatio](double accumulator, PairDI& alphaStemp){return accumulator + alphaStemp.first*PotentiationDepressionRatio;});
//     // if (alphaStimmulus+alphaBasal<0.0){
//     //     alphaStimmulus= (-alphaBasal);
//     // }
// }

bool ResourceSynapseSpine::ApplyAllTemporaryEffectsOnDepression(double expDecayFactor)
{
    //This is because depression updates instantly, so there is no need to source the decayHashMap for the factor as it will be identical for all changes 
    alphaStimmulus=std::accumulate(depressionAlphaTempAndCount.begin(), depressionAlphaTempAndCount.end(), alphaStimmulus, [expDecayFactor](double accumulator, double alphaStemp){return accumulator - alphaStemp*expDecayFactor;});
    //This is to avoid transition into negative weights through LTD. Also, for a single time-step the weight will be zero.
    if (alphaStimmulus+alphaBasal<0.0){
        alphaStimmulus= (-alphaBasal);
    }
    //The bool return is to keep track of plasticity events (for now not very useful, as still no distinction by branch)
    return (depressionAlphaTempAndCount.size()!=0);
}

void ResourceSynapseSpine::TickStimmulusCounts()
{
    //Increases the time count in the pairs of the potentiation vector (sort of like the trace in a sense)
    std::for_each(potentiationAlphaTempAndCount.begin(), potentiationAlphaTempAndCount.end(), [](PairDI& element) {element.second++;});
}

void ResourceSynapseSpine::CullStimmulusVectors()
{
    // for (int i=stimmulusEffectCount.size()-1; i>0; i--){//If exception, i>1//This may need to be changed later to use a vector for the indexes (or delete everything at the left of the last 0 count)
    //     if (stimmulusEffectCount.at(i)>=maxCount){
    //         stimmulusEffectCount.erase(std::next(stimmulusEffectCount.begin(), i));
    //         kStimmulusTempEffect.erase(std::next(kStimmulusTempEffect.begin(), i));
    //         nStimmulusTempEffect.erase(std::next(nStimmulusTempEffect.begin(), i));
    //     }
    // }

    for (std::list<PairDI>::reverse_iterator reverseIterator = potentiationAlphaTempAndCount.rbegin(); reverseIterator != potentiationAlphaTempAndCount.rend(); reverseIterator++)
    {
        //Here we iterate over the alpha temporary effects on reverse. That way the first time we find a time-obsolete effect, we delete everything before it.
        //This is because everything before the last (first in reverse) count that reached max is older, so they must have reached the max count too
        if (reverseIterator->second>=maxCount){
            potentiationAlphaTempAndCount.erase(potentiationAlphaTempAndCount.begin(), reverseIterator.base());
            break;
        }
    }
    //Single-timestep nature of this vector means we can just .clear()
    depressionAlphaTempAndCount.clear();
}
