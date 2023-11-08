//
// Created by Saif Ahmed on 27.06.21.
//
//
// Refactored by Antoni Bertolin on 14.06.23
//
#include "MonoDendriteSTDP.hpp"
#include <string>
#include <numeric>



MonoDendriteSTDP::MonoDendriteSTDP(GlobalSimInfo *infoGlobal) : Morphology(infoGlobal), stepWeights(false)
{
}

void MonoDendriteSTDP::Advect() {

    this->WeightDecay();

    this->ThetaDecay();
    signed long totalSpines{static_cast<signed long>(this->spineDataCoop.size())};
    //Recording plasticity events here makes no sense, because every time-step all thetas change and all weights change, and the framework is not event-based.

    // update cooperativity between spiker and spiker pairs -> avoids double counting when combined with loop that follows
    for (signed long spineID1 : std::ranges::views::iota(0, static_cast<signed long>(spikedSpinesId.size()))) {
        for (signed long spineID2 : std::ranges::views::iota(spineID1+1, static_cast<signed long>(spikedSpinesId.size()))) { //Not permutation loop, it is triangular. 
            this->UpdateCooperativity(this->spikedSpinesId.at(spineID1), this->spikedSpinesId.at(spineID2));
        }
    }

    // update cooperativity between spiker and non-spiker pairs
    for (const int spikedSpineId : this->spikedSpinesId) {
        for (signed long spineID : std::ranges::views::iota(0,totalSpines)) {
            if (!this->spikedSpines.at(spineID) && spikedSpineId != spineID) {//
                this->UpdateCooperativity(spikedSpineId, spineID);
            }
        }
    }

    // update for post-pre effects for spiker synapses
    for (const auto spikedSpineId: this->spikedSpinesId) {
        if (this->lastPostSpikeTime > 0 && this->integratePostSpike.at(spikedSpineId)) {
            this->UpdateLTD(spikedSpineId);
            this->integratePostSpike.at(spikedSpineId) = false;
        }
    }

    // update for pre-post effects for all synapses
    if (this->postSpiked){
        for (const CoopSpinePtr& spine: this->spineDataCoop) {
            if (spine->GetLastSpike() > 0 && this->integratePreSpike.at(spine->idInMorpho)) {
                this->UpdateLTP(spine->idInMorpho);
                this->integratePreSpike.at(spine->idInMorpho) = false;
            }
        }
        this->postSpiked=false;
    }

    this->Reset();
}
void MonoDendriteSTDP::NormalizeWeights() {
    if (this->weightNormalization == HardNormalization) {
        this->HardNormalize();
    } else if (this->weightNormalization == SoftMaxNormalization) {
        this->SoftMaxNormalize();
    }
}
void MonoDendriteSTDP::WeightDecay() {
    if (this->decayWeights) {
        std::for_each(baseSpineData.begin(), baseSpineData.end(), [this](BaseSpinePtr spine){
            spine->weight=(spine->GetWeightUncoupled() * weightExpDecay);
        });
    }
}
void MonoDendriteSTDP::ThetaDecay() {
    std::for_each(spineDataCoop.begin(), spineDataCoop.end(), [this](const CoopSpinePtr& spine){
        spine->SetTheta(spine->GetTheta() * thetaExpDecay);
    });
}

void MonoDendriteSTDP::RecordPostSpike() {
    Morphology::RecordPostSpike();
    this->lastPostSpikeTime = this->infoGlobal->dtTimestep * static_cast<double> (this->infoGlobal->timeStep);
    // STDP Analysis
    //this->postSpikes.push_back(this->lastPostSpikeTime);
    std::fill(this->integratePostSpike.begin(), this->integratePostSpike.end(), true);
}

void MonoDendriteSTDP::RecordExcitatoryPreSpike(int spikedSpineId) {
//This function is NOT DELAY COMPATIBLE (careful with the delays in synapse objects)
    Morphology::RecordExcitatoryPreSpike(spikedSpineId);
    this->spikedSpines.at(spikedSpineId) = true;//This does not seem to be correctly implemented
    this->spikedSpinesId.push_back(spikedSpineId);
    this->spineDataCoop.at(spikedSpineId)->SetLastSpike(static_cast<double>(this->infoGlobal->timeStep) * this->infoGlobal->dtTimestep); //only coop
    this->integratePreSpike.at(spikedSpineId) = true;
}

void MonoDendriteSTDP::SaveParameters(std::ofstream& wParameterFile, std::string neuronIdentificator) const {
    Morphology::SaveParameters(wParameterFile, neuronIdentificator);

    wParameterFile << neuronIdentificator<<"weightNormalization\t\t\t";
    if (this->weightNormalization == HardNormalization){
        wParameterFile<<"HardNormalization\n";
    }
    else if (this->weightNormalization == SoftMaxNormalization){
        wParameterFile<<"SoftMaxNormalization\n";
    }
    else if (this->weightNormalization == NOPNormalization){
        wParameterFile<<"NOPNormalization\n";
    }
    
    wParameterFile << neuronIdentificator<<"dendriticLength\t\t\t"<<std::to_string(this->dendriticLength);
    wParameterFile << "\t"<<"#Length of the dendritic arm (μm).\n";

    wParameterFile << neuronIdentificator<<"synapticGap\t\t\t"<<std::to_string(this->synapticGap);
    wParameterFile << "\t"<<"#Forced distance between synapses (μm).\n";

    wParameterFile << neuronIdentificator<<"heterosynapticThetaDecay\t\t"<<std::to_string(this->tauTheta);
    wParameterFile << "\t"<<"#Or tauTheta, decay constant of heterosynaptic effects in spines.\n";

    wParameterFile << neuronIdentificator<<"intersynapseDistanceDecay\t\t"<<std::to_string(this->lambdaDist);
    wParameterFile << "\t"<<"#Or lambdaDist, spatial decay constant of heterosynaptic boost between synapses.\n";

    wParameterFile << neuronIdentificator<<"intersynapseSpikeTimingDecay\t"<<std::to_string(this->tauDelay);
    wParameterFile << "\t"<<"#Or tauDelay, decay constant of heterosynaptic effects over inter-synapse spike timing difference.\n";

    wParameterFile << neuronIdentificator<<"preFactorLtp\t\t\t"<<std::to_string(this->preFactorLTP);
    wParameterFile << "\t"<<"#Base factor that is multiplied by the spatio-temporal effects in LTP. If set to zero, LTP will be zero. \"A\" equivalent\n";

    wParameterFile << neuronIdentificator<<"preFactorLtd\t\t\t"<<std::to_string(this->preFactorLTD);
    wParameterFile << "\t"<<"#Base factor that is multiplied by the spatio-temporal effects in LTD. If set to zero, LTD will be zero. \"A\" equivalent\n";

    wParameterFile << neuronIdentificator<<"incrLtp\t\t\t\t"<<std::to_string(this->incrementLTP);
    wParameterFile << "\t"<<"#Max possible increase in LTP due to cooperativity . \"I\" equivalent\n";

    wParameterFile << neuronIdentificator<<"decrLtd\t\t\t\t"<<std::to_string(this->decrementLTD);
    wParameterFile << "\t"<<"#Max possible LTD due to cooperativity.  \"D\" equivalent\n";

    wParameterFile << neuronIdentificator<<"baseLtp\t\t\t\t"<<std::to_string(this->baseLTP);
    wParameterFile << "\t"<<"#Default increase in weight per LTP check.  \"B\" equivalent\n";

    wParameterFile << neuronIdentificator<<"baseLtd\t\t\t\t"<<std::to_string(this->baseLTD);
    wParameterFile << "\t"<<"#Default decrease of weight per LTD check. \"B\" equivalent\n";
    
    wParameterFile << neuronIdentificator<<"weightDecay\t\t\t"<<std::boolalpha<<this->decayWeights<<std::noboolalpha<<"\t"<<std::setprecision(3)<<std::to_string(this->WeightDecayConstant)<<std::setprecision(6);
    wParameterFile<<"\t"<<"#The first bool activates the weight decay per timestep. The second number is the time constant on an exponential in seconds [exp(-dt/ctt)].\n";

    wParameterFile << neuronIdentificator<<"minmaxWeights\t\t\t"<<std::setprecision(3)<<std::to_string(this->minWeight)<<" "<<std::to_string(this->maxWeight)<<std::setprecision(6);
    wParameterFile<<"\t"<<"#Only relevant for HardNormalization and distributeWeights, the first number is the minimum weight in normalization, the second one the hard cap for weight.\n";
    wParameterFile << neuronIdentificator<<"alphaLtp\t\t\t\t"<<std::to_string(this->alpha);
    wParameterFile << "\t"<<"#Cooperativity decay for LTP.\n";
    wParameterFile << neuronIdentificator<<"betaLtd\t\t\t\t"<<std::to_string(this->beta);
    wParameterFile << "\t"<<"#Cooperativity decay for LTD.\n";

}

void MonoDendriteSTDP::LoadParameters(const std::vector<FileEntry>& morphologyParameters) {
    Morphology::LoadParameters(morphologyParameters);

    bool branchLengthInitiallized {false},
         synapticGapInitialized {false},
         tauThetaInitialized {false},
         lambdaDistInitialized {false},
         tauDelayInitialized {false},
         baseLtpInitialized {false},
         incrLtpInitialized {false},
         baseLtdInitialized {false},
         decrLtdInitialized {false},
         alphaInitialized{false},
         betaInitialized{false},
         normalizationFound {false};


    this->preFactorLTP = 1.0;
    this->preFactorLTD = 1.0;

    for (auto& [parameterName, parameterValues] : morphologyParameters) {

        if(parameterName.find("dendriticLength") != std::string::npos){
            this->dendriticLength = std::stod(parameterValues.at(0));
            branchLengthInitiallized = true;
        } else if (parameterName.find("synapticGap") != std::string::npos) {
            this->synapticGap = std::stod(parameterValues.at(0));
            synapticGapInitialized = true;
        } else if (parameterName.find("heterosynapticThetaDecay") != std::string::npos) {
            this->tauTheta = std::stod(parameterValues.at(0));
            tauThetaInitialized = true;
        } else if (parameterName.find("intersynapseDistanceDecay") != std::string::npos){
            this->lambdaDist = std::stod(parameterValues.at(0));
            lambdaDistInitialized = true;
        } else if (parameterName.find("intersynapseSpikeTimingDecay") != std::string::npos) {
            this->tauDelay = std::stod(parameterValues.at(0));
            tauDelayInitialized = true;
        } else if (parameterName.find("weightNormalization") != std::string::npos) {
            if (parameterValues.at(0) == IDstringNOPNormalization) {
                weightNormalization = NOPNormalization;
                normalizationFound = true;
            }
            else if (parameterValues.at(0) == IDstringHardNormalization) {
                weightNormalization = HardNormalization;
                normalizationFound = true;
            }
            else if (parameterValues.at(0) == IDstringSoftMaxNormalization) {
                weightNormalization = SoftMaxNormalization;
                normalizationFound = true;
            }
         }else if (parameterName.find("distributeWeights") != std::string::npos) {
            //This whole part is experimental, it seems it was not completely tested
            //As such, this is deprecated from publication
            if (parameterValues.at(0) == "true") {
                distributeWeights = true;
            } else if (parameterValues.at(0) == "step") {
                stepWeights = true;
                int cIdx = 1;
                try {
                    int weightSteps = std::stoi(parameterValues.at(cIdx++));
                    for (int k = 0; k < weightSteps; k++) {
                        this->weightStepBoundary.push_back(std::stol(parameterValues.at(cIdx++)));
                        this->weightStepValue.push_back(std::stod(parameterValues.at(cIdx++)));
                    }
                    this->currWightStepId = 0;
                } catch (...) {
                    throw "Issues with step weights..";
                }
            } else {
                try {
                    this->initialWeights = std::stod(parameterValues.at(1));
                } catch (...) {
                    throw "No initial weight read from file";
                }
            }
        } else if (parameterName.find("preFactorLtp") != std::string::npos) {
            this->preFactorLTP = std::stod(parameterValues.at(0));
        } else if (parameterName.find("preFactorLtd") != std::string::npos) {
            this->preFactorLTD = std::stod(parameterValues.at(0));
        } else if (parameterName.find("baseLtp") != std::string::npos) {
            this->baseLTP = std::stod(parameterValues.at(0));
            baseLtpInitialized = true;
        } else if (parameterName.find("baseLtd") != std::string::npos) {
            this->baseLTD = std::stod(parameterValues.at(0));
            baseLtdInitialized = true;
        } else if (parameterName.find("incrLtp") != std::string::npos) {
            this->incrementLTP = std::stod(parameterValues.at(0));
            incrLtpInitialized = true;
        } else if (parameterName.find("decrLtd") != std::string::npos) {
            this->decrementLTD = std::stod(parameterValues.at(0));
            decrLtdInitialized = true;
        }else if (parameterName.find("weightDecay") != std::string::npos) {
            this->decayWeights = {parameterValues.at(0)=="true"};
            this->WeightDecayConstant = std::stod(parameterValues.at(1));
            this->weightExpDecay=exp(-this->infoGlobal->dtTimestep/this->WeightDecayConstant);
        } else if (parameterName.find("alphaLtp") != std::string::npos) {
            this->alpha = std::stod(parameterValues.at(0));
            alphaInitialized = true;
        } else if (parameterName.find("betaLtd") != std::string::npos) {
            this->beta = std::stod(parameterValues.at(0));
            betaInitialized = true;
        } else if (parameterName.find("minmaxWeights") != std::string::npos) {
            this->minWeight = std::stod(parameterValues.at(0));
            this->maxWeight = std::stod(parameterValues.at(1));
        }

    }
    if (!branchLengthInitiallized){
        throw "Using heterosynaptic synapses without specifying dendriticLength is not allowed.";
    } else if (!synapticGapInitialized){
        throw "Using heterosynaptic synapses without specifying synapticGap is not allowed.";
    } else if (!tauThetaInitialized){
        throw "Using heterosynaptic synapses without specifying tau_theta is not allowed.";
    } else if (!lambdaDistInitialized){
        throw "Using heterosynaptic synapses without specifying lambdaDist is not allowed.";
    } else if (!tauDelayInitialized){
        throw "Using heterosynaptic synapses without specifying tau_delay is not allowed.";
    } else if (!baseLtdInitialized || !baseLtpInitialized || !incrLtpInitialized || !decrLtdInitialized){
        throw "some of the incr/decr ltp/ltd params not set";
    } else if (!alphaInitialized){
        throw "Using heterosynaptic synapses without specifying alpha is not allowed.";
    } else if (!betaInitialized){
        throw "Using heterosynaptic synapses without specifying beta is not allowed.";
    } else if (!normalizationFound){
        throw "No normalization found";
    }
//    this->posLo = this->synapticGap;
//    this->posHi = this->dendriticLength - this->synapticGap;
//    this->allocateDistal = false;
    this->nextPos =  this->synapticGap;
    this->spineIdGenerator = 0;
    this->thetaExpDecay=exp(-this->infoGlobal->dtTimestep/this->tauTheta);
}

void MonoDendriteSTDP::CheckParameters(const std::vector<FileEntry> &parameters) {
    Morphology::CheckParameters(parameters);
    for (auto& [parameterName, parameterValues] : parameters)  {

        if(parameterName.find("dendriticLength") != std::string::npos){
            if(!(this->dendriticLength == std::stod(parameterValues.at(0)))){
                throw "dendriticLength was not consistent in plasticity parameters";
            }
        } else if (parameterName.find("synapticGap") != std::string::npos) {
            if(!(this->synapticGap == std::stod(parameterValues.at(0)))){
                throw "synapticGap was not consistent in plasticity parameters";
            }
        } else if (parameterName.find("heterosynapticThetaDecay") != std::string::npos) {
            if(!(this->tauTheta == std::stod(parameterValues.at(0)))){
                throw "heterosynapticThetaDecay was not consistent in plasticity parameters";
            }
        } else if (parameterName.find("intersynapseDistanceDecay") != std::string::npos){
            if(!(this->lambdaDist == std::stod(parameterValues.at(0)))){
                throw "intersynapseDistanceDecay was not consistent in plasticity parameters";
            }
        } else if (parameterName.find("intersynapseSpikeTimingDecay") != std::string::npos) {
            if(!(this->tauDelay == std::stod(parameterValues.at(0)))){
                throw "intersynapseSpikeTimingDecay was not consistent in plasticity parameters";
            }
        } else if (parameterName.find("weightNormalization") != std::string::npos) {
            if (parameterValues.at(0) == IDstringNOPNormalization) {
                if(weightNormalization != NOPNormalization){
                    throw "weightNormalization was not consistent in plasticity parameters";
                }
            }
            else if (parameterValues.at(0) == IDstringHardNormalization) {
                if(weightNormalization != HardNormalization){
                    throw "weightNormalization was not consistent in plasticity parameters";
                }
            }
            else if (parameterValues.at(0) == IDstringSoftMaxNormalization) {
                if(weightNormalization != SoftMaxNormalization){
                    throw "weightNormalization was not consistent in plasticity parameters";
                }
            }
         }else if (parameterName.find("distributeWeights") != std::string::npos) {
            //This whole part is experimental, it seems it was not completely tested
            //As such, this is deprecated from publication
            if (parameterValues.at(0) == "true") {
                if (!distributeWeights){
                    throw "distributeWeights was not consistent in plasticity parameters";
                }
            } else if (parameterValues.at(0) == "step") {
                if (!stepWeights){
                    throw "distributeWeights was not consistent in plasticity parameters";
                }
            } else {
                if (initialWeights != std::stod(parameterValues.at(1))){
                    throw "distributeWeights was not consistent in plasticity parameters";
                }
            }
        } else if (parameterName.find("preFactorLtp") != std::string::npos) {
            if(!(this->preFactorLTP == std::stod(parameterValues.at(0)))){
                throw "preFactorLtp was not consistent in plasticity parameters";
            }
        } else if (parameterName.find("preFactorLtd") != std::string::npos) {
            if(!(this->preFactorLTD == std::stod(parameterValues.at(0)))){
                throw "preFactorLtd was not consistent in plasticity parameters";
            }
        } else if (parameterName.find("baseLtp") != std::string::npos) {
            if(!(this->baseLTP == std::stod(parameterValues.at(0)))){
                throw "baseLtp was not consistent in plasticity parameters";
            }
        } else if (parameterName.find("baseLtd") != std::string::npos) {
            if(!(this->baseLTD == std::stod(parameterValues.at(0)))){
                throw "baseLtd was not consistent in plasticity parameters";
            }
        } else if (parameterName.find("incrLtp") != std::string::npos) {
            if(!(this->incrementLTP == std::stod(parameterValues.at(0)))){
                throw "incrLtp was not consistent in plasticity parameters";
            }
        } else if (parameterName.find("decrLtd") != std::string::npos) {
            if(!(this->decrementLTD == std::stod(parameterValues.at(0)))){
                throw "decrLtd was not consistent in plasticity parameters";
            }
        }else if (parameterName.find("weightDecay") != std::string::npos) {
            if(this->decayWeights != (parameterValues.at(0)=="true")){
                throw "weightDecay was not consistent in plasticity parameters";
            } else if(!(this->WeightDecayConstant == std::stod(parameterValues.at(1)))){
                throw "weightDecay was not consistent in plasticity parameters";
            }
        } else if (parameterName.find("alphaLtp") != std::string::npos) {
            if(!(this->alpha == std::stod(parameterValues.at(0)))){
                throw "alphaLtp was not consistent in plasticity parameters";
            }
        } else if (parameterName.find("betaLtd") != std::string::npos) {
            if(!(this->beta == std::stod(parameterValues.at(0)))){
                throw "betaLtd was not consistent in plasticity parameters";
            }
        } else if (parameterName.find("minmaxWeights") != std::string::npos) {
            if(!(this->minWeight == std::stod(parameterValues.at(0)))){
                throw "minmaxWeights was not consistent in plasticity parameters";
            } else if(!(this->maxWeight == std::stod(parameterValues.at(1)))){
                throw "minmaxWeights was not consistent in plasticity parameters";
            }
        }

    }
}

void MonoDendriteSTDP::Reset() {
    this->NormalizeWeights();
    //std::fill(this->spikedSynapses.begin(),this->spikedSynapses.end(), false);
    std::fill(this->spikedSpines.begin(),this->spikedSpines.end(), false);
    this->spikedSpinesId.clear();
}
void MonoDendriteSTDP::HardNormalize() {
    for (BaseSpinePtr spine: this->baseSpineData) {
        spine->weight=(std::max(minWeight, std::min(maxWeight, spine->GetWeightUncoupled())));
    }
}

void MonoDendriteSTDP::SoftMaxNormalize() {

    //Softmax normalization (NNs version)
    // maxWeights = std::numeric_limits<double>::min();
    // for (auto& syn : this->baseSpineData) {
    //     maxWeights = std::max(maxWeights, syn->GetWeight());
    // }

    double weightSum = std::accumulate(this->baseSpineData.begin(), this->baseSpineData.end(), 0.0, [] (double weightSum, BaseSpinePtr synapse){
        return weightSum += std::exp(synapse->GetWeightUncoupled()); 
        });
    // for (auto& syn : this->baseSpineData) {
    //     sumWeights += std::exp(syn->GetWeight() - maxWeights);
    // }

    //It is not clear if the following lines are correct in Softmax normalization. There was no reference previously, so Toni assumed the normalization is the one done in NNs.
    //As for the multiplication by two, this is because the weights in Saif models are normally distributed between 0 and 2. This can be changed in a model with an extra entry in LP
    std::for_each(this->baseSpineData.begin(), this->baseSpineData.end(), [weightSum, this](BaseSpinePtr synapse){
        synapse->weight=((std::exp(synapse->GetWeightUncoupled())*this->softMaxMultiplier)/weightSum);
        });
}
BaseSpinePtr MonoDendriteSTDP::AllocateNewSynapse(BranchTargeting& branchTargeting) {

    std::uniform_real_distribution<double> distribution(0.0,2.0);

    // CoopSpinePtr newSpine = std::make_unique<CoopSynapseSpine>();
    this->spineDataCoop.push_back(std::make_unique<CoopSynapseSpine>());
    CoopSynapseSpine* newSpine = spineDataCoop.back().get();
    if (this->nextPos <= this->dendriticLength) { //This is very badly thought. The synapse allocation should throw, and not allow the programme to keep going with an unresolved issue
//        if (this->allocateDistal) {
//            newSynapse->distToSoma = this->posHi;
//            this->posHi -= this->synapticGap;
//        } else {
//            newSynapse->distToSoma = this->posLo;
//            this->posLo += this->synapticGap;
//        }
//
//        this->allocateDistal = !(this->allocateDistal;

        newSpine->SetDistToSoma( this->nextPos);
        newSpine->SetLastSpike(-200.0); // large negative value indicates no spikes of synapse during simulation
        newSpine->SetTheta(0);
        if (stepWeights) {
            if (static_cast<int>(spineDataCoop.size()) > weightStepBoundary.at(currWightStepId)) {
                currWightStepId++;
            }
            newSpine->weight=(weightStepValue.at(currWightStepId));
        } else {
            if (distributeWeights) {
                newSpine->weight=(distribution(generator));
            } else {
                newSpine->weight=(this->initialWeights); // assuming a range of weight between 0 and 2, weight is initialized to midpoint: 1
            }
        }
        //this->weightsSum += newSynapse->GetWeight();
        newSpine->idInMorpho=(this->spineIdGenerator++);
        // newSynapse->postNeuronId = ? // set in the Synapse object that calls for a new synapse
        // newSynapse->preNeuronId = ? // set in the Synapse object that calls for a new synapse
        this->baseSpineData.push_back(static_cast<BaseSpinePtr>(newSpine));

        this->nextPos += this->synapticGap;

        this->spikedSpines.push_back(false);
        this->integratePostSpike.push_back(false);
        this->integratePreSpike.push_back(false);
    } else {
        throw "MonoDendriteSTDP";
    }

    return static_cast<BaseSpinePtr>(newSpine);
}

void MonoDendriteSTDP::UpdateCooperativity(signed long spikerId, signed long neighborId) {
    CoopSynapseSpine* spiker = this->spineDataCoop.at(spikerId).get();
    CoopSynapseSpine* neighbor = this->spineDataCoop.at(neighborId).get();

    double hEffects = getDistanceEffects(spiker, neighbor) * getTimingEffects(spiker, neighbor) *spiker->GetWeightUncoupled() * neighbor->GetWeightUncoupled();
    spiker->AddToTheta(hEffects);
    neighbor->AddToTheta(hEffects);

//    if (hEffects != 0.0) {//OPTIMIZATION. Problems with heap allocation of the data (causes an OOM error). Needs to be moved to a file if necessary
////        pseudoCoop( spikerId, neighborId);
//
//        this->theta_changes.emplace_back(spikerId, hEffects);
//        this->theta_changes.emplace_back(neighborId, hEffects);
//    }
}
double MonoDendriteSTDP::aLTP(double theta) const {
    return baseLTP + incrementLTP * (1 - exp(-this->alpha * theta));
}
double MonoDendriteSTDP::aLTD(double theta) const {
    return (baseLTD - decrementLTD * (1 - exp(-this->beta * theta)));
}
double MonoDendriteSTDP::getDistanceEffects(const CoopSynapseSpine *const spineA, const CoopSynapseSpine *const spineB) const
{
    if (spineA == spineB) {
        return 0.0;
    } else {
        return exp(-std::abs(spineA->GetDistToSoma() - spineB->GetDistToSoma()) / this->lambdaDist);
    }
}
double MonoDendriteSTDP::getTimingEffects(const CoopSynapseSpine *const spineA, const CoopSynapseSpine *const spineB) const {
    if ((spineA == spineB) || (spineA->GetLastSpike() < 0 || spineB->GetLastSpike() < 0)) {
        return 0.0;
    } else {
        return exp(-std::abs(spineA->GetLastSpike() - spineB->GetLastSpike()) / this->tauDelay);
    }
}
// void MonoDendriteSTDP::pseudoCoop(signed long synId, signed long neighborId) {
//     CoopSynapseSpine* spiker = this->spineDataCoop.at(synId).get();
//     CoopSynapseSpine* neighbor = this->spineDataCoop.at(neighborId).get();
//     std::cout << "id 1: " << synId << ", id 2: " << neighborId << std::endl;
//     std::cout << "dist: " << abs(spiker->GetDistToSoma() - neighbor->GetDistToSoma()) << std::endl;
//     std::cout << "time: " << abs(spiker->GetLastSpike() - neighbor->GetLastSpike()) << std::endl;
//     std::cout << "dist effect: " << getDistanceEffects(spiker, neighbor);
//     std::cout << "time effect: " << getTimingEffects(spiker, neighbor);

// }