//
// Created by Antoni Bertolin on 14.06.23
//
#include "ResourceSynapseSpine.hpp"



std::vector<double> ResourceSynapseSpine::GetIndividualSynapticProfile() const
{
    std::vector<double> dataArray(5);
    dataArray.at(0) = this->branchId;
    dataArray.at(1) = this->weight;
    dataArray.at(2) = this->alphaResources;
    dataArray.at(3) = this->distanceFromNode;
    dataArray.at(4) = this->prePopId;
    return dataArray;
}

std::string ResourceSynapseSpine::GetIndividualSynapticProfileHeaderInfo() const{
    return std::string("{<branch ID>, <weight>, <alpha resources>, <distance to node>, <presynaptic population ID>}");
}

void ResourceSynapseSpine::ComputeAlphaResources(){
    //alphaResources=(kBasal+kStimulus)/(nBasal+nStimulus);
    alphaResources=alphaBasal+alphaStim;
    DecayAlphaResources();
}

void ResourceSynapseSpine::DecayAlphaResources(){
    // kStimulus*=kStimulusExpDecay;
    // nStimulus*=nStimulusExpDecay;
    alphaStim *= alphaStimulusExpDecay;
}

void ResourceSynapseSpine::ComputeWeight(double resourceFactor){
    this->weight=this->alphaResources*resourceFactor;
}

// void ResourceSynapseSpine::AddTempResourcesToSpine(double alphaStimulusInput){
//     // kStimulusTempEffect.push_back(kStimulusInput);
//     // nStimulusTempEffect.push_back(nStimulusInput);
//     //stimulusEffectCount.push_back(0);
//     depAlphaTempAndCount.push_back(alphaStimulusInput);
//     potAlphaTempCount.emplace_back(std::pair<double, int>(std::abs(alphaStimulusInput), 0));
// }

// bool ResourceSynapseSpine::ApplyAllTemporaryEffectsOnPostspike(const std::unordered_map<int, double>& STDPdecayMap) {
//     // kStimulus=std::accumulate(kStimulusTempEffect.begin(), kStimulusTempEffect.end(), kStimulus, [PotentiationDepressionRatio](double accumulator, double kStemp){return accumulator + kStemp*PotentiationDepressionRatio;});
//     // nStimulus=std::accumulate(nStimulusTempEffect.begin(), nStimulusTempEffect.end(), nStimulus, [PotentiationDepressionRatio](double accumulator, double nStemp){return accumulator + nStemp*PotentiationDepressionRatio;});
//     // kStimulusTempEffect.clear();
//     // nStimulusTempEffect.clear();
//     //stimulusEffectCount.clear();
//     //PotentiationDepressionRatio is to boost potentiation compared to depression, but also to input the current decay
//     //The next line goes over all the pairs in the potentiation vector and adds their stored alpha stimulus (stored there by pairing events) and adds them to the alphaStimulus variable
//     //In addition, it multiplies with the decay to postspike of those effects and multiplies by the potentiation/depression ratio (a parameter from .txt)
//     alphaStim=std::accumulate(potAlphaTempCount.begin(), potAlphaTempCount.end(), alphaStim, [this, STDPdecayMap](double accumulator, const std::pair<double, int>& alphaStemp){
//         return accumulator + alphaStemp.first**STDPdecayMap.at(alphaStemp.second);
//         });
//     // if (alphaStimulus+alphaBasal<0.0){
//     //     alphaStimulus= (-alphaBasal);
//     // }
//     //The bool return is to keep track of plasticity events (for now not very useful, as still no distinction by branch)
//     return (!potAlphaTempCount.empty());
// }

// // void ResourceSynapseSpine::ApplyAllTempEffectsOnConflictPotentiation(double PotentiationDepressionRatio)
// // {
// //     //CUrrently I have no way of properly calculating the decay of potentiation proper by calculating the area under the graph.
// //     alphaStimulus=std::accumulate(potentiationAlphaTempAndCount.begin(), potentiationAlphaTempAndCount.end(), alphaStimulus, [PotentiationDepressionRatio](double accumulator, std::pair<double, int>& alphaStemp){return accumulator + alphaStemp.first*PotentiationDepressionRatio;});
// //     // if (alphaStimulus+alphaBasal<0.0){
// //     //     alphaStimulus= (-alphaBasal);
// //     // }
// // }

// bool ResourceSynapseSpine::ApplyAllTemporaryEffectsOnDepression(double expDecayFactor) {
//     //This is because depression updates instantly, so there is no need to source the decayHashMap for the factor as it will be identical for all changes 
//     alphaStim=std::accumulate(depAlphaTempAndCount.begin(), depAlphaTempAndCount.end(), alphaStim, [expDecayFactor](double accumulator, double alphaStemp){
//         return accumulator - alphaStemp*expDecayFactor;
//         });
//     //This is to avoid transition into negative weights through LTD. Also, for a single time-step the weight will be zero.
//     alphaStim=std::min(alphaStim, -alphaBasal);
//     //The bool return is to keep track of plasticity events (for now not very useful, as still no distinction by branch)
//     return (!depAlphaTempAndCount.empty());
// }

// void ResourceSynapseSpine::TickStimulusCounts() {
//     //Increases the time count in the pairs of the potentiation vector (sort of like the trace in a sense)
//     std::for_each(potAlphaTempCount.begin(), potAlphaTempCount.end(), [](std::pair<double, int>& element) {
//         element.second++;
//     });
// }

// void ResourceSynapseSpine::CullStimulusVectors() {
//     // for (int i=stimulusEffectCount.size()-1; i>0; i--){//If exception, i>1//This may need to be changed later to use a vector for the indexes (or delete everything at the left of the last 0 count)
//     //     if (stimulusEffectCount.at(i)>=maxCount){
//     //         stimulusEffectCount.erase(std::next(stimulusEffectCount.begin(), i));
//     //         kStimulusTempEffect.erase(std::next(kStimulusTempEffect.begin(), i));
//     //         nStimulusTempEffect.erase(std::next(nStimulusTempEffect.begin(), i));
//     //     }
//     // }

//     for (std::list<std::pair<double, int>>::reverse_iterator reverseIterator = potAlphaTempCount.rbegin(); reverseIterator != potAlphaTempCount.rend(); reverseIterator++){
//         //Here we iterate over the alpha temporary effects on reverse. That way the first time we find a time-obsolete effect, we delete everything before it.
//         //This is because everything before the last (first in reverse) count that reached max is older, so they must have reached the max count too
//         if (reverseIterator->second>=maxCount){
//             potAlphaTempCount.erase(potAlphaTempCount.begin(), reverseIterator.base());
//             break;
//         }
//     }
//     //Single-timestep nature of this vector means we can just .clear()
//     depAlphaTempAndCount.clear();
// }
void ResourceSynapseSpine::AddTraceToAlpha(double computedTrace) {
    alphaStim = std::max(alphaStim + alphaStimBump*computedTrace, -alphaBasal);
}