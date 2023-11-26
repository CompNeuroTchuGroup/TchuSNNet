//
// Created by Antoni Bertolin on 14.06.23
//
#include "./AlphaResourceHSTDP.hpp"
#include "AlphaResourceHSTDP.hpp"

AlphaResourceHSTDP::AlphaResourceHSTDP(GlobalSimInfo *infoGlobal) : BranchedMorphology(infoGlobal) {
}

void AlphaResourceHSTDP::LoadParameters(const std::vector<FileEntry> &morphologyParameters) {

    // BranchedMorphology::LoadParameters(morphologyParameters);

    for (auto &[parameterName, parameterValues] : morphologyParameters) {
        if (parameterName.find("basalAlpha") != std::string::npos) {
            this->alphaBasal = std::stod(parameterValues.at(0));
        } else if (parameterName.find("alphaTau") != std::string::npos) {
            this->alphaStimulusTau      = std::stod(parameterValues.at(0));
            this->alphaStimulusExpDecay = std::exp(-this->infoGlobal->dtTimestep / this->alphaStimulusTau);
        } else if (parameterName.find("baseAlphaIncrease") != std::string::npos) {
            this->baseAlphaStimBump = std::stod(parameterValues.at(0));
        } else if (parameterName.find("omegaOffset") != std::string::npos) {
            this->omegaOffset = std::stod(parameterValues.at(0));
        } else if (parameterName.find("betaResourcePool") != std::string::npos) {
            this->betaResourcePool = std::stod(parameterValues.at(0));
            // } else if (parameterName.find("kernel_spatial_length") != std::string::npos){
            //     this->kernelGapNumber = static_cast<int>(std::stod(parameterValues.at(0))/this->synapticGap);
            // } else if (parameterName.find("kernel_temporal_length") != std::string::npos){
            //     this->timeKernelLength =
            //     static_cast<int>(std::stod(parameterValues.at(0))/this->infoGlobal->dtTimestep);
        } else if (parameterName.find("coopTau") != std::string::npos) {
            this->cooperativityTau = std::stod(parameterValues.at(0));
            this->coopExpDecay     = std::exp(-infoGlobal->dtTimestep / cooperativityTau);
        } else if (parameterName.find("coopProfile") != std::string::npos) {
            this->spaceProfileLambda = std::stod(parameterValues.at(0));
        } else if (parameterName.find("tauSTDP") != std::string::npos) {
            this->tauSTDP      = std::stod(parameterValues.at(0));
            this->STDPExpDecay = std::exp(-infoGlobal->dtTimestep / tauSTDP);
        } else if (parameterName.find("biasLTD") != std::string::npos) {
            this->biasLTD = std::stod(parameterValues.at(0));
            // } else if (parameterName.find("STDP_time_window") != std::string::npos){
            //     this->MaxCountSTDP = static_cast<int>(std::stod(parameterValues.at(0))/this->infoGlobal->dtTimestep);
            //     this->STDPDepressionCount=MaxCountSTDP;
        }
    }
    BranchedMorphology::LoadParameters(morphologyParameters); // Branchings are set up inside this call
    // Here branchings are already set up
    this->betaResourcePool /= branches.size();
    SetUpHashTable();
}

void AlphaResourceHSTDP::CheckParameters(const std::vector<FileEntry> &parameters) {
    BranchedMorphology::CheckParameters(parameters);
    for (auto &[parameterName, parameterValues] : parameters) {
        if (parameterName.find("basalAlpha") != std::string::npos) {
            if (!(this->alphaBasal == std::stod(parameterValues.at(0)))) {
                throw "basalAlpha was not consistent in plasticity model parameters.";
            }
        } else if (parameterName.find("alphaTau") != std::string::npos) {
            if (!(this->alphaStimulusTau == std::stod(parameterValues.at(0)))) {
                throw "alphaTau was not consistent in plasticity model parameters.";
            }
        } else if (parameterName.find("baseAlphaIncrease") != std::string::npos) {
            if (!(this->baseAlphaStimBump == std::stod(parameterValues.at(0)))) {
                throw "baseAlphaIncrease was not consistent in plasticity model parameters.";
            }
        } else if (parameterName.find("omegaOffset") != std::string::npos) {
            if (!(this->omegaOffset == std::stod(parameterValues.at(0)))) {
                throw "baseAlphaIncrease was not consistent in plasticity model parameters.";
            }
        } else if (parameterName.find("betaResourcePool") != std::string::npos) {
            if (!(this->betaResourcePool == std::stod(parameterValues.at(0)))) {
                throw "betaResourcePool was not consistent in plasticity model parameters.";
            }
            // } else if (parameterName.find("kernel_spatial_length") != std::string::npos){
            //     this->kernelGapNumber = static_cast<int>(std::stod(parameterValues.at(0))/this->synapticGap);
            // } else if (parameterName.find("kernel_temporal_length") != std::string::npos){
            //     this->timeKernelLength =
            //     static_cast<int>(std::stod(parameterValues.at(0))/this->infoGlobal->dtTimestep);
        } else if (parameterName.find("coopTau") != std::string::npos) {
            if (!(this->cooperativityTau == std::stod(parameterValues.at(0)))) {
                throw "coopTau was not consistent in plasticity model parameters.";
            }
        } else if (parameterName.find("coopProfile") != std::string::npos) {
            if (!(this->spaceProfileLambda == std::stod(parameterValues.at(0)))) {
                throw "coopProfile was not consistent in plasticity model parameters.";
            }
        } else if (parameterName.find("tauSTDP") != std::string::npos) {
            if (!(this->tauSTDP == std::stod(parameterValues.at(0)))) {
                throw "tauSTDP was not consistent in plasticity model parameters.";
            }
        } else if (parameterName.find("biasLTD") != std::string::npos) {
            if (!(this->biasLTD == std::stod(parameterValues.at(0)))) {
                throw "biasLTD was not consistent in plasticity model parameters.";
            }
            // } else if (parameterName.find("STDP_time_window") != std::string::npos){
            //     this->MaxCountSTDP = static_cast<int>(std::stod(parameterValues.at(0))/this->infoGlobal->dtTimestep);
            //     this->STDPDepressionCount=MaxCountSTDP;
        }
    }
}

void AlphaResourceHSTDP::SaveParameters(std::ofstream &wParameterFile, std::string neuronIdentificator) const {
    BranchedMorphology::SaveParameters(wParameterFile, neuronIdentificator);
    wParameterFile << neuronIdentificator << "basalAlpha\t\t\t" << std::to_string(this->alphaBasal);
    wParameterFile << "\t"
                   << "#Alpha at rest, where alpha decays towards\n";

    wParameterFile << neuronIdentificator << "alphaTau\t\t\t" << std::to_string(this->alphaStimulusTau);
    wParameterFile << " #secs\t"
                   << "#Decay constant of Alpha_stimulus\n";

    wParameterFile << neuronIdentificator << "baseAlphaIncrease\t\t" << std::to_string(this->baseAlphaStimBump);
    wParameterFile << "\t"
                   << "#Default alphaStimulus increase before applying spatial and temporal decays\n";

    wParameterFile << neuronIdentificator << "omegaOffset\t\t\t" << std::to_string(this->omegaOffset);
    wParameterFile << "\t"
                   << "#Offset factor in the weight definition.\n";

    wParameterFile << neuronIdentificator << "betaResourcePool\t\t"
                   << std::to_string(this->betaResourcePool * noBranches); //different var is for debugging branch generation
    wParameterFile << "\t"
                   << "#Multiplication factor of the definition of weight, representing the available 'total "
                      "resources'. Evenly split among branches\n";

    // wParameterFile <<
    // neuronIdentificator<<"kernel_spatial_length\t"<<std::to_string(this->kernelGapNumber*this->synapticGap);
    // wParameterFile << " #μm\t"<<"#Limit distance between two spines to be considered for synaptic spine pairing.\n";

    // wParameterFile <<
    // neuronIdentificator<<"kernel_temporal_length\t"<<std::to_string(this->timeKernelLength*this->infoGlobal->dtTimestep);//CHANGE
    // wParameterFile << " #secs\t"<<"#Maximum time length between two spikes in contiguous spines to be considered for
    // synaptic spine pairing\n";

    wParameterFile << neuronIdentificator << "coopTau\t\t\t" << std::to_string(this->cooperativityTau); // CHANGE
    wParameterFile << " #secs\t"
                   << "#Time decay constant for the alpha stimulus increase in the pairing kernel\n";

    wParameterFile << neuronIdentificator << "coopProfile\t\t\t" << std::to_string(this->spaceProfileLambda); // CHANGE
    wParameterFile << " #μm\t"
                   << "#Exponential decay that characterizes the spatial profile of cooperativity\n";

    wParameterFile << neuronIdentificator << "tauSTDP\t\t\t" << std::to_string(this->tauSTDP); // CHANGE
    wParameterFile << " #secs\t"
                   << "#Exponential decay constant for the STDP kernel. Affects both pre and post synaptic traces\n";

    // wParameterFile <<
    // neuronIdentificator<<"STDP_time_window\t\t"<<std::to_string(this->MaxCountSTDP*this->infoGlobal->dtTimestep);//CHANGE
    // wParameterFile << " #secs\t"<<"#Max time where STDP potentiation/depression can happen\n";

    wParameterFile << neuronIdentificator << "biasLTD\t\t\t" << std::to_string(this->biasLTD); // CHANGE
    wParameterFile << "\t"
                   << "#Factor that biases the resting STDP towards potentiation or depression. 1 is symmetric, lower "
                      "is potentiating, higher is depressing\n";

    wParameterFile << "##### The weight of this model is defined as wI=beta*alphaI/(omega+branch-sum(alpha)), where "
                      "alphaI= alphaBasal + alphaStimulus*exp(-dt/alphaStimTau) \n";
}

int AlphaResourceHSTDP::CreateBranch(std::vector<int> anteriorBranches) {
    int branchId{this->GenerateBranchId()};
    if (anteriorBranches.empty()) {
        this->alphaBranches.emplace_back(AlphaBranch(anteriorBranches,this->synapticGap, this->branchLength,  branchId, STDPExpDecay,
                        coopExpDecay,alphaStimulusExpDecay,baseAlphaStimBump,alphaBasal, betaResourcePool, omegaOffset)); // This vector should be sorted by ID by default (tested).
        this->branches.push_back(static_cast<Branch *>(&this->alphaBranches.back()));
    } else {
        int branchId{this->GenerateBranchId()};
        this->alphaBranches.emplace_back(
            AlphaBranch(this->synapticGap, this->branchLength, branchId, STDPExpDecay,
                        coopExpDecay, alphaStimulusExpDecay, baseAlphaStimBump, alphaBasal, betaResourcePool, omegaOffset)); // This vector should be sorted by ID by default (tested).
        this->branches.push_back(static_cast<Branch *>(&this->alphaBranches.back()));
    }
    return branchId;
}

void AlphaResourceHSTDP::SetUpHashTable() {
    // Statement 2 of second for loop is done thinking about time steps being in the order of magnitude of 0.1 ms, while
    // the gap is 1 um. If we make the kernel triangular, every "step" is one gap in space and in time it can be any
    // amount of timesteps. To calculate the step relationship, we do:

    spatialProfile.resize(std::round(branchLength / synapticGap));
    for (int spaceIndex : std::ranges::views::iota(1, static_cast<int>(spatialProfile.size()))) {
        spatialProfile.at(spaceIndex) = std::exp((-synapticGap * spaceIndex) / spaceProfileLambda);
    }
}

void AlphaResourceHSTDP::Advect() {
    // std::for_each(rTBranches.begin(), rTBranches.end(), [this](const ResourceTraceBranch* const branch){
    //     UpdateCoopTrace(branch);
    // });
    if (this->postSpiked) {
        // If post spike, apply all stimms on positive mode (remember the coded function in spines) with the decay from
        // STDP pot count.
        // Use the count in the effects of synapses for the actual decay for STDP, but the branch vector for detecting
        // the updatable ones WITH DECAY (of alpha, STDP-like)
        std::for_each(alphaBranches.begin(), alphaBranches.end(),
                      [this](AlphaBranch &branch) { branch.ApplyTracesOnSpinesLTP(); });

    } else if (CheckIfPreSpikeHappened()) { // Checks if all are empty or some are not. MAy be redundant if frequency is
                                            // high enough.
        ApplyPreSpikePerturbations();
        std::for_each(alphaBranches.begin(), alphaBranches.end(), [this](AlphaBranch &branch) {
            branch.ApplyTracesOnSpinesLTD(this->GetPostSynapticTrace(), GetLTDBias());
        });
    }
    // ApplyEffects();
    Reset();
}


bool AlphaResourceHSTDP::CheckIfPreSpikeHappened() {
    return std::any_of(PAR_UNSEQ, alphaBranches.begin(), alphaBranches.end(),
                       [](const AlphaBranch &branch) { return !branch.spikedSpinesInTheBranch.empty(); });
}


void AlphaResourceHSTDP::ApplyCoopTraceSpatialProfile(int branchSpineID, AlphaBranch &const currentBranch) {
    // This function will evaluate all spines within kernel distance of the synapse spine and queue the alpha stimm
    // effects for each pairing found (queue waiting for postspike) All reference declarations are to reduce indexing
    // times in containers 

    std::vector<double>::iterator forwardBegin = currentBranch.cooperativityTraces.begin() + (branchSpineID + 1);
    std::vector<double>::iterator forwardEnd   = currentBranch.cooperativityTraces.end();
    std::transform(PAR_UNSEQ, forwardBegin, forwardEnd, spatialProfile.begin(), forwardBegin, std::plus<double>());
    
    std::vector<double>::reverse_iterator reverseBegin =
        currentBranch.cooperativityTraces.rbegin() + (currentBranch.cooperativityTraces.size() - branchSpineID);
    std::vector<double>::reverse_iterator reverseEnd = currentBranch.cooperativityTraces.rend();
    std::transform(PAR_UNSEQ, reverseBegin, reverseEnd, spatialProfile.begin(), reverseBegin, std::plus<double>());
}
void AlphaResourceHSTDP::ApplyPreSpikePerturbations() {
        // branch.preSynapticTraces.at(branchSpinePosition) += 1;
    // ApplyCoopTraceSpatialProfile(branchSpinePosition, branch);
    std::for_each(alphaBranches.begin(), alphaBranches.end(), [this](AlphaBranch& branch){
        for (int spinePosID: branch.spikedSpinesInTheBranch){
            branch.preSynapticTraces.at(spinePosID) += 1;
            ApplyCoopTraceSpatialProfile(spinePosID, branch);
        }
    });
}


void AlphaResourceHSTDP::Reset() {
    BranchedMorphology::Reset();
    // this->postSpiked=false;
    // Wrapper plus clearing some of the vectors. Last method to run in chronological order, where we call the ticks and
    // the general upkeep Check increase beta resources here? Not necessary anymore Clear both sets in every branch DONE
    // REVIEW if there are any remaining things to put in this function
    // ClearSynapseSets();
    std::for_each(alphaBranches.begin(), alphaBranches.end(), [this](AlphaBranch &branch) { branch.ComputeWeights(); });
    DecayAllTraces();
    // DeleteEffects();
}

void AlphaResourceHSTDP::DecayAllTraces() {
    std::for_each(alphaBranches.begin(), alphaBranches.end(), [](AlphaBranch &branch) { branch.DecayAllTraces(); });
    // std::for_each(resourceSpineData.begin(), resourceSpineData.end(), [](ResourceSpinePtr spine){
    //     spine->TickStimulusCounts();
    // });
    postSynapticTrace *= STDPExpDecay; // If it is unbound, no bool checks
    // if (STDPDepressionCount<MaxCountSTDP){
    //     STDPDepressionCount++;
    // }
}



void AlphaResourceHSTDP::RecordPostSpike() {
    this->totalPostSpikes++;
    this->postSpiked = true;
    postSynapticTrace += 1;
}

void AlphaResourceHSTDP::RecordExcitatoryPreSpike(BaseSpinePtr spinePtr) {
    // This function is NOT DELAY COMPATIBLE (careful with the delays in synapse objects)
    //POINTER CASTING IS FREE PERFORMANCE WISE
    // Here only record, afterwards we do the checks
    // Not going down the virtual path because inefficient
    // AlphaSynapseSpine *synapseSpine = alphaSpines.at(spikedSpineId);
    // AlphaBranch&          branch       = alphaBranches.at(synapseSpine->branchId);
    // int                   branchSpinePosition{synapseSpine->branchPositionId};
    alphaBranches.at(static_cast<AlphaSpinePtr>(spinePtr)->branchId).spikedSpinesInTheBranch.push_back(static_cast<AlphaSpinePtr>(spinePtr)->branchPositionId);
    // branch.preSynapticTraces.at(branchSpinePosition) += 1;
    // ApplyCoopTraceSpatialProfile(branchSpinePosition, branch);
    // branch->cooperativityTraces.at(synapseSpine->branchPositionId)+=1;
    this->totalPreSpikes++;
}

void AlphaResourceHSTDP::PostConnectSetUp() {
    BranchedMorphology::PostConnectSetUp();
    if (distributeWeights) {
        std::uniform_real_distribution<double> weightDistribution(-alphaBasal, alphaBasal);
        for (AlphaBranch& branch : alphaBranches){
            for (double& alphaStim: branch.spineAlphaStims){
                alphaStim=(weightDistribution(generator));
            }
        }
    }
}

BaseSpinePtr AlphaResourceHSTDP::AllocateNewSynapse(BranchTargeting &branchTarget) {
    // here I have to set the maxcount of the spine to maxCount too
    // Here sum over the branches.????
    // And cast the proper pointers to the proper baseSpineData vectors.
    int branch{AllocateBranch(branchTarget)};
    int position{PopSynapseSlotFromBranch(branchTarget)};
    alphaBranches.at(branch).alphaSpines.push_back(AlphaSynapseSpine());
    AlphaSynapseSpine *newSpine = &alphaBranches.at(branch).alphaSpines.back();
    this->alphaSpines.push_back(newSpine);
    // this->weightsSum += newSynapse->GetWeight();
    newSpine->idInMorpho=this->baseSpineData.size();//this->spineIdGenerator++
    // Branch
    newSpine->branchPositionId=position;
    newSpine->branchId=branch;

    branches.at(branch)->synapseSlotClosedIndex.push_back(position);
    // branches.at(branch)->morphoSynapseIDs.push_back(newSynapse->GetIdInMorpho());
    // branches.at(branch)->synapseSlotToMorphoIDMap[position]=newSynapse->GetIdInMorpho();
    // Storage (other)
    this->baseSpineData.push_back(static_cast<BaseSpinePtr>(newSpine));
    this->branchedSpineData.push_back(static_cast<BranchedSpinePtr>(newSpine));

    return this->baseSpineData.back();
}

std::vector<double> AlphaResourceHSTDP::GetOverallSynapticProfile() const {
    /*
     * returned array organised as follows:
     * item 1: average synaptic weight
     * item 2: totalLTD Events
     * item 3: totalLTP Events
     * item 4: average plasticity events
     * */
    std::vector<double> dataArray(4);
    size_t              noSpines{this->baseSpineData.size()};
    // size_t sizeOfSpineData {this->baseSpineData.size()};
    //  double weightSum = std::accumulate(this->baseSpineData.begin(), this->baseSpineData.end(), 0.0, [] (double
    //  accumulator, const BaseSpinePtr syn) {
    //                                      return accumulator + syn->GetWeight();
    //                                      });
    // CalcMorphoPlasticityEvents();

    dataArray.at(0) = std::accumulate(this->baseSpineData.begin(), this->baseSpineData.end(), 0.0,
                                      [](double accumulator, const BaseSpinePtr syn) {
                                          return accumulator + syn->GetWeightUncoupled();
                                      }) /
                      noSpines;
    dataArray.at(1) = this->totalPostSpikes;
    dataArray.at(2) = this->totalPreSpikes;
    dataArray.at(3) = noSpines;
    return dataArray;
}

std::string AlphaResourceHSTDP::GetOverallSynapticProfileHeaderInfo() const {
    return std::string("{<average weight>, <total post spikes> ,<total pre spikes>, <total number of spines>}");
}

// void BranchedResourceHeteroSTDP::CalcMorphoPlasticityEvents()
// {
//     totalLTDEvents = std::accumulate(this->branches.begin(), this->branches.end(), 0,
//                                        [] (int accumulator, const BranchPtr branch) { return accumulator +
//                                        branch->LTDevents; });
//     totalLTPEvents = std::accumulate(this->branches.begin(), this->branches.end(), 0,
//                                        [] (int accumulator, const BranchPtr branch) { return accumulator +
//                                        branch->LTPevents; });
//     totalPlasticityEvents = totalLTDEvents + totalLTPEvents;
// }
