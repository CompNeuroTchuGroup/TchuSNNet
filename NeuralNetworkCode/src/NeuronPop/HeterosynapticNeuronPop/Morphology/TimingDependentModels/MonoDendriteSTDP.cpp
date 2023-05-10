//
// Created by Saif Ahmed on 27.06.21.
//

#include "MonoDendriteSTDP.hpp"
#include <string>
#include <numeric>

MonoDendriteSTDP::MonoDendriteSTDP(GlobalSimInfo* info): Morphology(info), integratePostSpike(false), initialWeights(1.0), stepWeights(false) {
}

void MonoDendriteSTDP::advect() {

    this->WeightDecay();

    this->thetaDecay();

    //Recording plasticity events here makes no sense, because every time-step all thetas change and all weights change, and the framework is not event-based.

    // update cooperativity between spiker and spiker pairs -> avoids double counting when combined with loop that follows
    for (unsigned long i = 0; i < this->spikedSynapsesId.size(); ++i) {
        for (unsigned long j = i+1; j < this->spikedSynapsesId.size(); ++j) {
            this->updateCooperativity(this->spikedSynapsesId.at(i), this->spikedSynapsesId.at(j));
        }
    }

    // update cooperativity between spiker and non-spiker pairs
    for (const auto& spikedSynapseId : this->spikedSynapsesId) {
        for (unsigned long j = 0; j < this->synapseDataCoop.size(); ++j) {
            if (!this->spikedSynapses[j] && spikedSynapseId != j) {
                this->updateCooperativity(spikedSynapseId, j);
            }
        }
    }

    // update for post-pre effects for spiker synapses
    for (const auto& spikedSynapseId: this->spikedSynapsesId) {
        if (this->lastPostSpikeTime > 0 && this->integratePostSpike.at(spikedSynapseId)) {
            this->updateLTD(spikedSynapseId);
            this->integratePostSpike.at(spikedSynapseId) = false;

        }
    }

    // update for pre-post effects for all synapses
    if (this->postSpiked){
        for (const auto& syn: this->synapseDataCoop) {
            if (syn->GetLastSpike() > 0 && this->integratePreSpike.at(syn->GetIdInMorpho())) {
                this->updateLTP(syn->GetIdInMorpho());
                this->integratePreSpike.at(syn->GetIdInMorpho()) = false;
            }
        }
        this->postSpiked=false;
    }

    this->Reset();
}

void MonoDendriteSTDP::thetaDecay() {
    for (const std::shared_ptr<CoopSynapseSpine>& syn: this->synapseDataCoop) {
        syn->SetTheta(syn->GetTheta() * this->thetaExpDecay);
    }
}

void MonoDendriteSTDP::RecordPostSpike() {
    Morphology::RecordPostSpike();
    this->lastPostSpikeTime = this->info->dt * static_cast<double> (this->info->time_step);
    // STDP Analysis
    //this->postSpikes.push_back(this->lastPostSpikeTime);
    std::fill(this->integratePostSpike.begin(), this->integratePostSpike.end(), true);
}

void MonoDendriteSTDP::RecordExcitatoryPreSpike(int spikedSynapseId) {
//This function is NOT DELAY COMPATIBLE (careful with the delays in synapse objects)
    Morphology::RecordExcitatoryPreSpike(spikedSynapseId);
    this->spikedSynapses.at(spikedSynapseId) = true;//This does not seem to be correctly implemented
    this->spikedSynapsesId.push_back(spikedSynapseId);
    this->synapseDataCoop.at(spikedSynapseId)->SetLastSpike(static_cast<double> (this->info->time_step) * this->info->dt); //only coop
    this->integratePreSpike.at(spikedSynapseId) = true;
}

void MonoDendriteSTDP::SaveParameters(std::ofstream *stream, std::string neuronPreId) {
    Morphology::SaveParameters(stream, neuronPreId);

    *stream << neuronPreId<<"_morphology_weight_normalization\t";
    if (this->weightNormalization == HardNormalization){
        *stream<<"HardNormalization\n";
    }
    else if (this->weightNormalization == SoftMaxNormalization){
        *stream<<"SoftMaxNormalization\n";
    }
    else if (this->weightNormalization == NOPNormalization){
        *stream<<"NOPNormalization\n";
    }
    
    *stream << neuronPreId<<"_dendritic_length\t\t\t"<<std::to_string(this->dendriticLength);
    *stream << "\t"<<"#Length of the dendritic arm (μm).\n";

    *stream << neuronPreId<<"_synaptic_gap\t\t\t\t"<<std::to_string(this->synapticGap);
    *stream << "\t"<<"#Forced distance between synapses (μm).\n";

    *stream << neuronPreId<<"_heterosynaptic_theta_decay\t"<<std::to_string(this->tauTheta);
    *stream << "\t"<<"#Or tauTheta, decay constant of heterosynaptic effects in spines.\n";

    *stream << neuronPreId<<"_intersynapse_distance_decay\t"<<std::to_string(this->lambdaDist);
    *stream << "\t"<<"#Or lambdaDist, spatial decay constant of heterosynaptic boost between synapses.\n";

    *stream << neuronPreId<<"_intersynapse_spike_timing_decay\t"<<std::to_string(this->tauDelay);
    *stream << "\t"<<"#Or tauDelay, decay constant of heterosynaptic effects over inter-synapse spike timing difference.\n";

    *stream << neuronPreId<<"_pre_factor_ltp\t\t\t"<<std::to_string(this->preFactorLTP);
    *stream << "\t"<<"#Base factor that is multiplied by the spatio-temporal effects in LTP. If set to zero, LTP will be zero. \"A\" equivalent\n";

    *stream << neuronPreId<<"_pre_factor_ltd\t\t\t"<<std::to_string(this->preFactorLTD);
    *stream << "\t"<<"#Base factor that is multiplied by the spatio-temporal effects in LTD. If set to zero, LTD will be zero. \"A\" equivalent\n";

    *stream << neuronPreId<<"_incr_ltp\t\t\t\t"<<std::to_string(this->incr_ltp);
    *stream << "\t"<<"#Max possible increase in LTP due to cooperativity . \"I\" equivalent\n";

    *stream << neuronPreId<<"_decr_ltd\t\t\t\t"<<std::to_string(this->decr_ltd);
    *stream << "\t"<<"#Max possible LTD due to cooperativity.  \"D\" equivalent\n";

    *stream << neuronPreId<<"_base_ltp\t\t\t\t"<<std::to_string(this->base_ltp);
    *stream << "\t"<<"#Default increase in weight per LTP check.  \"B\" equivalent\n";

    *stream << neuronPreId<<"_base_ltd\t\t\t\t"<<std::to_string(this->base_ltd);
    *stream << "\t"<<"#Default decrease of weight per LTD check. \"B\" equivalent\n";
    
    *stream << neuronPreId<<"_morphology_weight_decay\t"<<std::boolalpha<<this->decayWeights<<std::noboolalpha<<"\t"<<std::to_string(this->WeightDecayConstant);
    *stream<<"\t"<<"#The first bool activates the weight decay per timestep. The second number is the time constant on an exponential in seconds [exp(-dt/ctt)].\n";

    *stream << neuronPreId<<"_morphology_min-max_weights\t"<<std::to_string(this->minWeight)<<"\t"<<std::to_string(this->maxWeight);
    *stream<<"\t"<<"#Only relevant for HardNormalization and distribute_weights, the first number is the minimum weight in normalization, the second one the hard cap for weight.\n";
}

void MonoDendriteSTDP::LoadParameters(std::vector<std::string> *input) {
    Morphology::LoadParameters(input);

    std::string name;
    std::vector<std::string> values;

    bool dendriteInitialized {false},
         synapticGapInitialized {false},
         tauThetaInitialized {false},
         lambdaDistInitialized {false},
         tauDelayInitialized {false},
         base_ltpInitialized {false},
         incr_ltpInitialized {false},
         base_ltdInitialized {false},
         decr_ltdInitialized {false},
         normalizationFound {false};


    this->preFactorLTP = 1.0;
    this->preFactorLTD = 1.0;

    for (auto & it : *input) {
        SplitString(&it,&name,&values);

        if(name.find("dendritic_length") != std::string::npos){
            this->dendriticLength = std::stod(values.at(0));
            dendriteInitialized = true;
        } else if (name.find("synaptic_gap") != std::string::npos) {
            this->synapticGap = std::stod(values.at(0));
            synapticGapInitialized = true;
        } else if (name.find("heterosynaptic_theta_decay") != std::string::npos) {
            this->tauTheta = std::stod(values.at(0));
            tauThetaInitialized = true;
        } else if (name.find("intersynapse_distance_decay") != std::string::npos){
            this->lambdaDist = std::stod(values.at(0));
            lambdaDistInitialized = true;
        } else if (name.find("intersynapse_spike_timing_decay") != std::string::npos) {
            this->tauDelay = std::stod(values.at(0));
            tauDelayInitialized = true;
        } else if (name.find("weight_normalization") != std::string::npos) {
            if (values.at(0) == str_NOPNormalization) {
                weightNormalization = NOPNormalization;
                normalizationFound = true;
            }
            else if (values.at(0) == str_HardNormalization) {
                weightNormalization = HardNormalization;
                normalizationFound = true;
            }
            else if (values.at(0) == str_SoftMaxNormalization) {
                weightNormalization = SoftMaxNormalization;
                normalizationFound = true;
            }
         }else if (name.find("distribute_weights") != std::string::npos) {
            //This whole part is experimental, it seems it was not completely tested
            //As such, this is deprecated from publication
            if (values.at(0) == "true") {
                distributeWeights = true;
            } else if (values.at(0) == "step") {
                stepWeights = true;
                int cIdx = 1;
                try {
                    int weightSteps = std::stoi(values.at(cIdx++));
                    for (int k = 0; k < weightSteps; k++) {
                        this->weightStepBoundary.push_back(std::stoul(values.at(cIdx++)));
                        this->weightStepValue.push_back(std::stod(values.at(cIdx++)));
                    }
                    this->currWightStepId = 0;
                } catch (...) {
                    bool x =false;
                    assertm(x == true, "Issues with step weights..");
                }
            } else {
                try {
                    this->initialWeights = std::stod(values.at(1));
                } catch (...) {
                }
            }
        } else if (name.find("pre_factor_ltp") != std::string::npos) {
            this->preFactorLTP = std::stod(values.at(0));
        } else if (name.find("pre_factor_ltd") != std::string::npos) {
            this->preFactorLTD = std::stod(values.at(0));
        } else if (name.find("base_ltp") != std::string::npos) {
            this->base_ltp = std::stod(values.at(0));
            base_ltpInitialized = true;
        } else if (name.find("base_ltd") != std::string::npos) {
            this->base_ltd = std::stod(values.at(0));
            base_ltdInitialized = true;
        } else if (name.find("incr_ltp") != std::string::npos) {
            this->incr_ltp = std::stod(values.at(0));
            incr_ltpInitialized = true;
        } else if (name.find("decr_ltd") != std::string::npos) {
            this->decr_ltd = std::stod(values.at(0));
            decr_ltdInitialized = true;
        }else if (name.find("weight_decay") != std::string::npos) {
            this->decayWeights = {values.at(0)=="true"};
            this->WeightDecayConstant = std::stod(values.at(1));
            this->weightExpDecay=exp(-this->info->dt/this->WeightDecayConstant);
        } else if (name.find("min-max_weights") != std::string::npos) {
            this->minWeight = std::stod(values.at(0));
            this->maxWeight = std::stod(values.at(1));
        }

    }

    assertm(dendriteInitialized == true, "Using heterosynaptic synapses without specifying dendritic_length is not allowed.");
    assertm(synapticGapInitialized == true, "Using heterosynaptic synapses without specifying synaptic_gap is not allowed.");

    // TODO: Update assert messages to appropriate
    assertm(tauThetaInitialized == true, "Using heterosynaptic synapses without specifying dendritic_length is not allowed.");
    assertm(lambdaDistInitialized == true, "Using heterosynaptic synapses without specifying synaptic_gap is not allowed.");
    assertm(tauDelayInitialized == true, "Using heterosynaptic synapses without specifying synaptic_gap is not allowed.");
    assertm(base_ltdInitialized == true && base_ltpInitialized == true && incr_ltpInitialized == true && decr_ltdInitialized == true,
            "some of the incr/decr ltp/ltd params not set");


//    this->posLo = this->synapticGap;
//    this->posHi = this->dendriticLength - this->synapticGap;
//    this->allocateDistal = false;
    this->nextPos =  this->synapticGap;
    this->synapseIdGenerator = 0;
    this->thetaExpDecay=exp(-this->info->dt/this->tauTheta);
}

void MonoDendriteSTDP::Reset()
{
    this->NormalizeWeights();
    //std::fill(this->spikedSynapses.begin(),this->spikedSynapses.end(), false);
    std::fill(this->spikedSynapses.begin(),this->spikedSynapses.end(), false);
    this->spikedSynapsesId.clear();
}

BaseSpinePtr MonoDendriteSTDP::AllocateNewSynapse(const HeteroCurrentSynapse& synapse) {

    std::uniform_real_distribution<double> distribution(0.0,2.0);

    std::shared_ptr<CoopSynapseSpine> newSynapse = std::make_shared<CoopSynapseSpine>();

    if (this->nextPos <= this->dendriticLength) { //This is very badly thought. The synapse allocation should throw, and not allow the programme to keep going with an unresolved issue
//        if (this->allocateDistal) {
//            newSynapse->distToSoma = this->posHi;
//            this->posHi -= this->synapticGap;
//        } else {
//            newSynapse->distToSoma = this->posLo;
//            this->posLo += this->synapticGap;
//        }
//
//        this->allocateDistal = !this->allocateDistal;

        newSynapse->SetDistToSoma( this->nextPos);
        newSynapse->SetLastSpike(-200.0); // large negative value indicates no spikes of synapse during simulation
        newSynapse->SetTheta(0);
        if (stepWeights) {
            if (synapseDataCoop.size() > weightStepBoundary.at(currWightStepId)) {
                currWightStepId++;
            }
            newSynapse->SetWeight(weightStepValue.at(currWightStepId));
        } else {
            if (distributeWeights) {
                std::random_device rd;
                std::mt19937 generator(rd()); // rd() provides a random seed
                newSynapse->SetWeight(distribution(generator));
            } else {
                newSynapse->SetWeight(this->initialWeights); // assuming a range of weight between 0 and 2, weight is initialized to midpoint: 1
            }
        }
        //this->weightsSum += newSynapse->GetWeight();
        newSynapse->SetIdInMorpho(this->synapseIdGenerator++);
        // newSynapse->postNeuronId = ? // set in the Synapse object that calls for a new synapse
        // newSynapse->preNeuronId = ? // set in the Synapse object that calls for a new synapse
        this->baseSynapseData.push_back(static_cast<BaseSpinePtr>(newSynapse));
        this->synapseDataCoop.push_back(newSynapse);
        this->nextPos += this->synapticGap;

        this->spikedSynapses.push_back(false);
        this->integratePostSpike.push_back(false);
        this->integratePreSpike.push_back(false);
    } else {
        throw;
    }

    return static_cast<BaseSpinePtr>(newSynapse);
}

void MonoDendriteSTDP::updateCooperativity(unsigned long spikerId, unsigned long neighborId) {
    CoopSynapseSpine* spiker = this->synapseDataCoop.at(spikerId).get();
    CoopSynapseSpine* neighbor = this->synapseDataCoop.at(neighborId).get();

    double hEffects = getDistanceEffects(spiker, neighbor);
    hEffects *= getTimingEffects(spiker, neighbor);
    hEffects *= spiker->GetWeight() * neighbor->GetWeight();
    spiker->AddToTheta(hEffects);
    neighbor->AddToTheta(hEffects);

//    if (hEffects != 0.0) {//OPTIMIZATION. Problems with heap allocation of the data (causes an OOM error). Needs to be moved to a file if necessary
////        pseudoCoop( spikerId, neighborId);
//
//        this->theta_changes.emplace_back(spikerId, hEffects);
//        this->theta_changes.emplace_back(neighborId, hEffects);
//    }
}

void MonoDendriteSTDP::pseudoCoop(unsigned long synId, unsigned long neighborId) {
    CoopSynapseSpine* spiker = this->synapseDataCoop.at(synId).get();
    CoopSynapseSpine* neighbor = this->synapseDataCoop.at(neighborId).get();
    std::cout << "id 1: " << synId << ", id 2: " << neighborId << std::endl;
    std::cout << "dist: " << abs(spiker->GetDistToSoma() - neighbor->GetDistToSoma()) << std::endl;
    std::cout << "time: " << abs(spiker->GetLastSpike() - neighbor->GetLastSpike()) << std::endl;
    std::cout << "dist effect: " << getDistanceEffects(spiker, neighbor);
    std::cout << "time effect: " << getTimingEffects(spiker, neighbor);

}