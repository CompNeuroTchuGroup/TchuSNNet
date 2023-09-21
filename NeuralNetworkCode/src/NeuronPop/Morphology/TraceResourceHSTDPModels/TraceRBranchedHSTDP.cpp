//
// Created by Antoni Bertolin on 14.06.23
//
#include "./TraceRBranchedHSTDP.hpp"
#include "TraceRBranchedHSTDP.hpp"

TraceRBranchedHSTDP::TraceRBranchedHSTDP(GlobalSimInfo* infoGlobal):BranchedMorphology(infoGlobal) {
}

void TraceRBranchedHSTDP::LoadParameters(const std::vector<FileEntry>& morphologyParameters) {
    
    // BranchedMorphology::LoadParameters(morphologyParameters);

    for (auto& [parameterName, parameterValues] : morphologyParameters) {
        if (parameterName.find("basalAlpha") != std::string::npos){
            this->alphaBasal = std::stod(parameterValues.at(0));
        } else if (parameterName.find("alphaTau") != std::string::npos){
            this->alphaStimulusTau = std::stod(parameterValues.at(0));
            this->alphaStimulusExpDecay = std::exp(-this->infoGlobal->dtTimestep/this->alphaStimulusTau);
        } else if (parameterName.find("baseAlphaIncrease") != std::string::npos){
            this->baseAlphaStimBump = std::stod(parameterValues.at(0));
        } else if (parameterName.find("omegaOffset") != std::string::npos){
            this->omegaOffset = std::stod(parameterValues.at(0));
        } else if(parameterName.find("betaResourcePool") != std::string::npos){
            this->betaResourcePool = std::stod(parameterValues.at(0));
        // } else if (parameterName.find("kernel_spatial_length") != std::string::npos){
        //     this->kernelGapNumber = static_cast<int>(std::stod(parameterValues.at(0))/this->synapticGap);
        // } else if (parameterName.find("kernel_temporal_length") != std::string::npos){
        //     this->timeKernelLength = static_cast<int>(std::stod(parameterValues.at(0))/this->infoGlobal->dtTimestep);
        } else if (parameterName.find("coopTau") != std::string::npos){
            this->cooperativityTau = std::stod(parameterValues.at(0));
            this->coopExpDecay =std::exp(-infoGlobal->dtTimestep/cooperativityTau);
        } else if (parameterName.find("coopProfile") != std::string::npos){
            this->spaceProfileLambda = std::stod(parameterValues.at(0));
        } else if (parameterName.find("tauSTDP") != std::string::npos){
            this->tauSTDP = std::stod(parameterValues.at(0));
            this->STDPExpDecay = std::exp(-infoGlobal->dtTimestep/tauSTDP);
        } else if (parameterName.find("biasLTD") != std::string::npos){
            this->biasLTD = std::stod(parameterValues.at(0));
        // } else if (parameterName.find("STDP_time_window") != std::string::npos){
        //     this->MaxCountSTDP = static_cast<int>(std::stod(parameterValues.at(0))/this->infoGlobal->dtTimestep);
        //     this->STDPDepressionCount=MaxCountSTDP;
        }
    }
    BranchedMorphology::LoadParameters(morphologyParameters);//Branchings are set up inside this call
    //Here branchings are already set up
    this->betaResourcePool/=branches.size();
    SetUpHashTable();
}

void TraceRBranchedHSTDP::CheckParameters(const std::vector<FileEntry> &parameters) {
    BranchedMorphology::CheckParameters(parameters);
    for (auto& [parameterName, parameterValues] : parameters) {
        if (parameterName.find("basalAlpha") != std::string::npos){
            if(!(this->alphaBasal == std::stod(parameterValues.at(0)))){
                throw "basalAlpha was not consistent in plasticity model parameters.";
            }
        } else if (parameterName.find("alphaTau") != std::string::npos){
            if(!(this->alphaStimulusTau == std::stod(parameterValues.at(0)))){
                throw "alphaTau was not consistent in plasticity model parameters.";
            }
        } else if (parameterName.find("baseAlphaIncrease") != std::string::npos){
            if(!(this->baseAlphaStimBump == std::stod(parameterValues.at(0)))){
                throw "baseAlphaIncrease was not consistent in plasticity model parameters.";
            }
        } else if (parameterName.find("omegaOffset") != std::string::npos){
            if(!(this->omegaOffset == std::stod(parameterValues.at(0)))){
                throw "baseAlphaIncrease was not consistent in plasticity model parameters.";
            }
        } else if(parameterName.find("betaResourcePool") != std::string::npos){
            if(!(this->betaResourcePool == std::stod(parameterValues.at(0)))){
                throw "betaResourcePool was not consistent in plasticity model parameters.";
            }
        // } else if (parameterName.find("kernel_spatial_length") != std::string::npos){
        //     this->kernelGapNumber = static_cast<int>(std::stod(parameterValues.at(0))/this->synapticGap);
        // } else if (parameterName.find("kernel_temporal_length") != std::string::npos){
        //     this->timeKernelLength = static_cast<int>(std::stod(parameterValues.at(0))/this->infoGlobal->dtTimestep);
        } else if (parameterName.find("coopTau") != std::string::npos){
            if(!(this->cooperativityTau == std::stod(parameterValues.at(0)))){
                throw "coopTau was not consistent in plasticity model parameters.";
            }
        } else if (parameterName.find("coopProfile") != std::string::npos){
            if(!(this->spaceProfileLambda == std::stod(parameterValues.at(0)))){
                throw "coopProfile was not consistent in plasticity model parameters.";
            }
        } else if (parameterName.find("tauSTDP") != std::string::npos){
            if(!(this->tauSTDP == std::stod(parameterValues.at(0)))){
                throw "tauSTDP was not consistent in plasticity model parameters.";
            }
        } else if (parameterName.find("biasLTD") != std::string::npos){
            if(!(this->biasLTD == std::stod(parameterValues.at(0)))){
                throw "biasLTD was not consistent in plasticity model parameters.";
            }
        // } else if (parameterName.find("STDP_time_window") != std::string::npos){
        //     this->MaxCountSTDP = static_cast<int>(std::stod(parameterValues.at(0))/this->infoGlobal->dtTimestep);
        //     this->STDPDepressionCount=MaxCountSTDP;
        }
    }
}

void TraceRBranchedHSTDP::SaveParameters(std::ofstream& wParameterFile, std::string neuronIdentificator) const {
    BranchedMorphology::SaveParameters(wParameterFile, neuronIdentificator);
    wParameterFile << neuronIdentificator<<"basalAlpha\t\t\t"<<std::to_string(this->alphaBasal);
    wParameterFile << "\t"<<"#Alpha at rest, where alpha decays towards\n";

    wParameterFile << neuronIdentificator<<"alphaTau\t\t\t"<<std::to_string(this->alphaStimulusTau);
    wParameterFile << " #secs\t"<<"#Decay constant of Alpha_stimulus\n";

    wParameterFile << neuronIdentificator<<"baseAlphaIncrease\t\t"<<std::to_string(this->baseAlphaStimBump);
    wParameterFile << "\t"<<"#Default alphaStimulus increase before applying spatial and temporal decays\n";

    wParameterFile << neuronIdentificator<<"omegaOffset\t\t\t"<<std::to_string(this->omegaOffset);
    wParameterFile << "\t"<<"#Offset factor in the weight definition.\n";

    wParameterFile << neuronIdentificator<<"betaResourcePool\t\t"<<std::to_string(this->betaResourcePool);
    wParameterFile << "\t"<<"#Multiplication factor of the definition of weight, representing the available 'total resources'. Evenly split among branches\n";

    // wParameterFile << neuronIdentificator<<"kernel_spatial_length\t"<<std::to_string(this->kernelGapNumber*this->synapticGap);
    // wParameterFile << " #μm\t"<<"#Limit distance between two spines to be considered for synaptic spine pairing.\n";

    // wParameterFile << neuronIdentificator<<"kernel_temporal_length\t"<<std::to_string(this->timeKernelLength*this->infoGlobal->dtTimestep);//CHANGE
    // wParameterFile << " #secs\t"<<"#Maximum time length between two spikes in contiguous spines to be considered for synaptic spine pairing\n";

    wParameterFile << neuronIdentificator<<"coopTau\t\t\t"<<std::to_string(this->cooperativityTau);//CHANGE
    wParameterFile << " #secs\t"<<"#Time decay constant for the alpha stimulus increase in the pairing kernel\n";

    wParameterFile << neuronIdentificator<<"coopProfile\t\t\t"<<std::to_string(this->spaceProfileLambda);//CHANGE
    wParameterFile << " #μm\t"<<"#Exponential decay that characterizes the spatial profile of cooperativity\n";

    wParameterFile << neuronIdentificator<<"tauSTDP\t\t\t"<<std::to_string(this->tauSTDP);//CHANGE
    wParameterFile << " #secs\t"<<"#Exponential decay constant for the STDP kernel. Affects both pre and post synaptic traces\n";
    
    // wParameterFile << neuronIdentificator<<"STDP_time_window\t\t"<<std::to_string(this->MaxCountSTDP*this->infoGlobal->dtTimestep);//CHANGE
    // wParameterFile << " #secs\t"<<"#Max time where STDP potentiation/depression can happen\n";

    wParameterFile << neuronIdentificator<<"biasLTD\t\t\t"<<std::to_string(this->biasLTD);//CHANGE
    wParameterFile << "\t"<<"#Factor that biases the resting STDP towards potentiation or depression. 1 is symmetric, lower is potentiating, higher is depressing\n";

    wParameterFile <<"##### The weight of this model is defined as wI=beta*alphaI/(omega+branch-sum(alpha)), where alphaI= alphaBasal + alphaStimulus*exp(-dt/alphaStimTau) \n";
}

void TraceRBranchedHSTDP::SetUpBranchings(int remainingBranchingEvents, std::vector<int> anteriorBranches) {
    //This is a recursive function that sets up the branched dendritic tree and is generalized for 0 branchings (1 branch). This function has been unit tested by Antoni.
    remainingBranchingEvents-=1;
    //First call is done with an empty int vector
    for (int index : std::ranges::views::iota(0,2)) {
        (void)index;
        int branchId{this->GenerateBranchId()};
        this->rTBranches.emplace_back(std::make_unique<ResourceTraceBranch>(ResourceTraceBranch(this->synapticGap, this->branchLength, anteriorBranches, branchId, STDPExpDecay, coopExpDecay)));//This vector should be sorted by ID by default (tested).
        this->branches.push_back(static_cast<Branch*>(this->rTBranches.back().get()));
        //Constructor here
        if(remainingBranchingEvents>0){
            std::vector<int> anteriorBranchesCopy(anteriorBranches);
            anteriorBranchesCopy.push_back(branchId);
            this->SetUpBranchings(remainingBranchingEvents, anteriorBranchesCopy);
        }
    }
}

void TraceRBranchedHSTDP::SetUpHashTable() {
//Statement 2 of second for loop is done thinking about time steps being in the order of magnitude of 0.1 ms, while the gap is 1 um. 
//If we make the kernel triangular, every "step" is one gap in space and in time it can be any amount of timesteps. To calculate the step relationship, we do:
    // if ((static_cast<double>(timeKernelLength)/kernelGapNumber)<1.0){
    //     throw "The timeKernel lenght is assumed to be longer than the spatialKernel length."; //The timeKernel lenght is always thought to be longer than the spatial one. 1 ms timing for 10 um is very one-sided, and that is equilibrium
    // }
    // timeStepsPerSynGap = static_cast<double>(timeKernelLength)/kernelGapNumber;
    // kernelRadiusArea=static_cast<double>(timeKernelLength)+timeStepsPerSynGap;
    // //The two loops generate a sort of ladder in the kernel (geometrically)
    // for (int spaceIndex : std::ranges::views::iota(0, kernelGapNumber+1)){
    //     std::unordered_map<int, double> tempHashMap{};
    //     for (int timeIndex : std::ranges::views::iota (0, timeKernelLength-(timeStepsPerSynGap*(spaceIndex-1))+1)){ 
    //         tempHashMap[timeIndex]=std::exp(-(synapticGap*spaceIndex)/spaceProfileLambda)*std::exp(-(infoGlobal->dtTimestep*timeIndex)/cooperativityTau);
    //         //Always access time, space later (always will be more times than gaps)
    //     }
    //     kernelHashTable[spaceIndex] = tempHashMap;
    // }

    // for (int STDPindex : std::ranges::views::iota (0, MaxCountSTDP+1)){
    //     DecayHashTableSTDP[STDPindex]=std::exp(-(STDPindex*infoGlobal->dtTimestep)/tauSTDP);
    // }
    spatialProfile.resize(std::round(branchLength/synapticGap));
    for (int spaceIndex : std::ranges::views::iota(1, static_cast<int>(spatialProfile.size()))){
        spatialProfile.at(spaceIndex) = std::exp((-synapticGap*spaceIndex)/spaceProfileLambda);
    }
}

void TraceRBranchedHSTDP::Advect() {
    // std::for_each(rTBranches.begin(), rTBranches.end(), [this](const ResourceTraceBranch* const branch){
    //     UpdateCoopTrace(branch);
    // });
    if (this->postSpiked){
         //If post spike, apply all stimms on positive mode (remember the coded function in spines) with the decay from STDP pot count. 
        //Use the count in the effects of synapses for the actual decay for STDP, but the branch vector for detecting the updatable ones
        //WITH DECAY (of alpha, STDP-like)
        std::for_each(rTBranches.begin(), rTBranches.end(), [this](RTBranchPtr& branch){branch->ApplyTracesOnSpinesLTP();});
        // for (const ResourceTraceBranch* const branch : rTBranches){
        //     for (ResourceSpinePtr spine : branch->rBranchSpineData){
        //         if (spine==nullptr){
        //             continue;
        //         } else {
        //             if(spine->ApplyAllTemporaryEffectsOnPostspike(DecayHashTableSTDP)){
        //             }
        //         }
        //     }
        // }
    } else if (CheckIfPreSpikeHappened()){ //Checks if all are empty or some are not. MAy be redundant if frequency is high enough.
        std::for_each(rTBranches.begin(), rTBranches.end(), [this](RTBranchPtr& branch){branch->ApplyTracesOnSpinesLTD(this->GetPostSynapticTrace(), this->GetLTDBias());});
    }
    // ApplyEffects();
    Reset();
}


// void TraceRBranchedHSTDP::ApplyEffects() { //Called after pairings
//     //If first precount is bigger than twice the STDP depression count, depression. Otherwise potentiation.
    
// }

// void TraceRBranchedHSTDP::UpdateCoopTrace(const ResourceTraceBranch* const branch) {//These are the spiked neurons on synapseDataIndexes 
//     //Here we call SetTempEffects if __it__ happens
//     std::for_each(branch->spikedSpinesInTheBranch.begin(), branch->spikedSpinesInTheBranch.end(), [branch, this](int spineIDinBranch){
//         ApplyCoopTraceSpatialProfile(spineIDinBranch, branch);
//     });
// }

bool TraceRBranchedHSTDP::CheckIfPreSpikeHappened() {
    return std::any_of(PAR_UNSEQ,rTBranches.begin(), rTBranches.end(), [](const RTBranchPtr& branch){return !branch->spikedSpinesInTheBranch.empty();});
}

// bool BranchedResourceHeteroSTDP::CheckIfThereIsPairing(RBranchPtr branch, int synapseIDinBranch)
// {
//     //This function raises an overflow flag in compilation due to operations before casting to ptrdiff_t (long long). There should never be overflow from a design standpoint
//     return std::count_if(std::max(branch->triggerCount.begin(), std::next(branch->triggerCount.begin(),synapseIDinBranch-kernelGapNumber)),std::min(branch->triggerCount.end(), std::next(branch->triggerCount.begin(),synapseIDinBranch+kernelGapNumber+1)), [this](int pairingCounter){return pairingCounter<this->timeKernelLength;})>2;
// }

void TraceRBranchedHSTDP::ApplyCoopTraceSpatialProfile(int branchSpineID, ResourceTraceBranch* const currentBranch) {
    //This function will evaluate all spines within kernel distance of the synapse spine and queue the alpha stimm effects for each pairing found (queue waiting for postspike)
    //All reference declarations are to reduce indexing times in containers
    //int spinePositionBranchIndex;//,absDistance,timeStepDifference;
    // double alphaStimulusEffect;
    //std::set<int>& kernelizedSynapses = currentBranch->updatedSynapseSpines;
    //int branchSlots{static_cast<int>(currentBranch->branchSlots)};
    // ResourceSynapseSpine* spikedSpine = currentBranch->rBranchSpineData.at(branchSpineID);
    std::vector<double>::iterator forwardBegin = currentBranch->cooperativityTraces.begin() + (branchSpineID+1);
    std::vector<double>::iterator forwardEnd = currentBranch->cooperativityTraces.end();
    std::transform(PAR_UNSEQ, forwardEnd, spatialProfile.begin(), forwardBegin, std::plus<double>());
    //Foward loop
    // for (int positionIndex : std::ranges::views::iota(1,branchSlots-branchSpineID)){//We loop from contiguous spine to the end of the branch, as represented by branchID-branchID+branchSlots-1=branchSlots-1, last index
    //     spinePositionBranchIndex = branchSpineID+positionIndex;
    //     // alphaStimulusEffect = baseAlphaStimBump;
    //     currentBranch->cooperativityTraces.at(spinePositionBranchIndex)+=spatialProfile.at(positionIndex);
    //     // if (spinePositionBranchIndex<0){//Condition to avoid illegal indexing
    //     //     continue;
    //     // } else if (spinePositionBranchIndex >= branchSlots) {//Condition to avoid illegal indexing
    //     //     break;
    //     // } else if (spinePositionBranchIndex == branchSpineID){
    //     //     spikedSpine->AddTempResourcesToSpine(alphaStimulusEffect);
    //     //     updatedSpines.push_back(branchSpineID);
    //     //     continue;//This essentially does regular STDP
    //     // } else {
    //     //     ResourceSpinePtr neighbourSynapseSpine = currentBranch->rBranchSpineData.at(spinePositionBranchIndex);
    //     //     if (neighbourSynapseSpine == nullptr) {
    //     //         continue; //To jump empty synapse slots
    //     //     }
    //     //     absDistance = std::abs(positionIndex);
    //     //     timeStepDifference=currentBranch->triggerCount.at(spinePositionBranchIndex);
    //     //     //timeStepDifference == triggerCount of first spine
    //     //     // if (timeStepDifference<=STDPDepressionCount){//Indicates depression region with no conflict. If the spike is far away, but too far for depression, nothing will happen. If a postspike happens, these flags are ignored.
    //     //     //     //The boolean is <= to classify the pairing where the postspike and first prespike are on the same timestep as depression.
    //     //     //     // centralSynapseSpine->SetDepressionFlag(true);
    //     //     //     // neighbourSynapseSpine->SetDepressionFlag(true);
    //     //     //     //alphaStimulusEffect *= -DecayHashTableSTDP.at(STDPDepressionCount-timeStepDifference);//This is to apply a decay equivalent to the STDP kernel to the first pre-spike with depression (*-1) (calcium accumulation oriented)
    //     //     // } else { //Knowing depression region, indicates conflict that is skewed towards potentiation according to STDP kernel
    //     //     //     if (timeStepDifference<=2*STDPDepressionCount){//We set the boolean to equal too to solve the conflict in favor of potentiation, as depression has the other fringe case
    //     //     //         centralSynapseSpine->SetPotentiationFlag(true);
    //     //     //         neighbourSynapseSpine->SetPotentiationFlag(true);
    //     //     //         alphaStimulusEffect *= DecayHashTableSTDP.at(timeStepDifference-STDPDepressionCount);//This is to apply a decay equivalent to the STDP kernel to the first pre-spike  (calcium accumulation oriented)
    //     //     //     } else {
    //     //     //         centralSynapseSpine->SetDepressionFlag(true);
    //     //     //         neighbourSynapseSpine->SetDepressionFlag(true);
    //     //     //     }//These flags are only here to solve 
    //     //     //}
    //     //     // if (STDPDepressionCount<MaxCountSTDP){
    //     //     //     alphaStimulusEffect=(-alphaStimulusEffect);
    //     //     // }
    //     //     if(timeStepDifference+timeStepsPerSynGap*absDistance<=kernelRadiusArea){ //&& (kernelizedSynapses.find(synapsePositionIndexInBranch) == kernelizedSynapses.end())){        //triggercount+distance gap*spaceTimeStepRelation<maxCountTime + 1*spaceTimeStepRelation (what should constrain the triangular matrix)
    //     //         //The postspiked condition is because the STDP kernel is discretized, in the case the postspike happens at the same time as the pairing of a synapse, the pairing does not happen (as the potentiation and depression would counteract themselves)
    //     //         //Avoiding the check of whether the neighbour was already updated assumes that no two spines will fire the same dt.
    //     //         //This assumption is necessary due to implementation limitations.
    //     //         alphaStimulusEffect *= CallKernelHashTable(absDistance, timeStepDifference);
    //     //         spikedSpine->AddTempResourcesToSpine(alphaStimulusEffect); //Removing this line (plus updated alpha effects) will change the model from pairing to cooperativity only.
    //     //         neighbourSynapseSpine->AddTempResourcesToSpine(alphaStimulusEffect);
    //     //         updatedSpines.push_back(spinePositionBranchIndex);
    //     //         updatedSpines.push_back(branchSpineID);
    //     //     }
    //     // }
    // }
    std::vector<double>::reverse_iterator reverseBegin = currentBranch->cooperativityTraces.rbegin() + (currentBranch->cooperativityTraces.size()-branchSpineID);
    std::vector<double>::reverse_iterator reverseEnd = currentBranch->cooperativityTraces.rend();
    std::transform(PAR_UNSEQ,reverseBegin, reverseEnd, spatialProfile.begin(), reverseBegin, std::plus<double>());
    // for (int positionIndex : std::ranges::views::iota(1,branchSpineID+1)){ //And then use branchSynapseID+gapindex
    // //We loop from contiguous spine until spine 0, from the index being branchID-branchID
    //     spinePositionBranchIndex = branchSpineID-positionIndex;
    //     // alphaStimulusEffect = baseAlphaStimBump;
    //     currentBranch->cooperativityTraces.at(spinePositionBranchIndex)+=spatialProfile.at(positionIndex);
    // }
    //Reverse loop
    //kernelizedSynapses.insert(branchSynapseID);//Used to avoid double potentiation in same timestep
    //Here for every synapse inside the synapse's kernel that has an active counter (count!=countMax) we get the time kernel and then apply the space kernel
    //Take unto account the synaptic GAP and DT!!! This should be done elsewhere
}
// double TraceRBranchedHSTDP::CallKernelHashTable(int distanceToCenterInGaps) {
//     return spatialProfile.at(distanceToCenterInGaps);
// }
// double TraceRBranchedHSTDP::CallKernelHashTable(int distanceToCenterInGaps, int timeDifference) {
//     //Matrix is symmetric!! Only have to do half!
//     //Here we tabulate for every distance from center a equivalency between counts and effects (reverse as counts go up to maxCount)
//     //The pair consists of ns and ks. The function has to obtain the count from synapse and distance?
//     //OR the hashtable contains the ns and ks multipliers (for now they are the same)
//     return kernelHashTable.at(distanceToCenterInGaps).at(timeDifference);
// }

void TraceRBranchedHSTDP::Reset() {
    BranchedMorphology::Reset();
    // this->postSpiked=false;
    //Wrapper plus clearing some of the vectors. Last method to run in chronological order, where we call the ticks and the general upkeep
    //Check increase beta resources here? Not necessary anymore
    //Clear both sets in every branch DONE
    //REVIEW if there are any remaining things to put in this function
    // ClearSynapseSets();
    std::for_each(rTBranches.begin(), rTBranches.end(), [this](RTBranchPtr& branch){
        ComputeWeights(branch.get());
    });
    DecayAllTraces();
    // DeleteEffects();
}

void TraceRBranchedHSTDP::ComputeAlphas(const ResourceTraceBranch* const branch) {
    //Here we just need to apply all delta alphas, decay them (before makes more sense, the bump has delta t zero.), sum the result to stationary alpha, then update alpha sums? Yes
    std::for_each(branch->rBranchSpineData.begin(), branch->rBranchSpineData.end(), [](ResourceSynapseSpine* const spine){
        if (spine != nullptr) {
            spine->ComputeAlphaResources();
        }
    });    
}

void TraceRBranchedHSTDP::ComputeWeights(ResourceTraceBranch* const branch) { //This one is the one we call for every branch {
    ComputeAlphaSums(branch);
    branch->resourceFactor=betaResourcePool/(omegaOffset+branch->alphaTotalSum);//IMPORTANT, if we make beta non-constant, beta must be referenced from the branch
    std::for_each(branch->rBranchSpineData.begin(), branch->rBranchSpineData.end(), [branch](ResourceSynapseSpine* const spine){
        if (spine != nullptr) {
            spine->ComputeWeight(branch->resourceFactor);
        }
    });
}

void TraceRBranchedHSTDP::ComputeAlphaSums(ResourceTraceBranch* const branch) {
    ComputeAlphas(branch);
    branch->alphaTotalSum=std::accumulate(branch->rBranchSpineData.begin(), branch->rBranchSpineData.end(), 0.0, [](double accumulator, const ResourceSynapseSpine* const spine) {
        if (spine != nullptr) {
            return accumulator + spine->GetAlphaResources();
        } else {
            return accumulator;
    }});
}
// void TraceRBranchedHSTDP::DeleteEffects() {
//     std::for_each(resourceSpineData.begin(), resourceSpineData.end(), [](ResourceSpinePtr spine){
//         spine->CullStimulusVectors();
//         spine->SetUpdatedFlag(false);
//     });
// }
void TraceRBranchedHSTDP::DecayAllTraces(){
    std::for_each(rTBranches.begin(), rTBranches.end(), [](RTBranchPtr& branch){
        branch->DecayAllTraces();
    });
    // std::for_each(resourceSpineData.begin(), resourceSpineData.end(), [](ResourceSpinePtr spine){
    //     spine->TickStimulusCounts();
    // });
    postSynapticTrace*=STDPExpDecay;//If it is unbound, no bool checks 
    // if (STDPDepressionCount<MaxCountSTDP){
    //     STDPDepressionCount++;
    // }
}

// void TraceRBranchedHSTDP::ClearSynapseSets() {
//     //Clearing the sets used for depression in this function (single timestep sets)
//     std::for_each(rTBranches.begin(), rTBranches.end(), [](const ResourceTraceBranch* const branch){
//         branch->updatedSpines.clear();
//     });
// }

void TraceRBranchedHSTDP::RecordPostSpike() {
    this->totalPostSpikes++;
    this->postSpiked = true;
    postSynapticTrace+=1;
}

void TraceRBranchedHSTDP::RecordExcitatoryPreSpike(int spikedSpineId) {
    //This function is NOT DELAY COMPATIBLE (careful with the delays in synapse objects)
    //Here only record, afterwards we do the checks
    //Not going down the virtual path because inefficient
    ResourceSynapseSpine* synapseSpine = resourceSpineData.at(spikedSpineId).get();
    ResourceTraceBranch* branch = rTBranches.at(synapseSpine->GetBranchId()).get();
    int branchSpinePosition{synapseSpine->GetBranchPositionId()};
    branch->spikedSpinesInTheBranch.push_back(branchSpinePosition);
    branch->preSynapticTraces.at(branchSpinePosition)+=1;
    ApplyCoopTraceSpatialProfile(branchSpinePosition, branch);
    // branch->cooperativityTraces.at(synapseSpine->GetBranchPositionId())+=1;
    this->totalPreSpikes++;
}

void TraceRBranchedHSTDP::PostConnectSetUp() {
    BranchedMorphology::PostConnectSetUp();
    if (distributeWeights){
        std::uniform_real_distribution<double> weightDistribution(-alphaBasal, alphaBasal);
        for(ResourceSpinePtr& spine : resourceSpineData){
            spine->SetAlphaStim(weightDistribution(generator));
        }
    }
}

BaseSpinePtr TraceRBranchedHSTDP::AllocateNewSynapse(const BranchTargeting& branchTarget) {
    //here I have to set the maxcount of the spine to maxCount too 
    //Here sum over the branches.????
    //And cast the proper pointers to the proper baseSpineData vectors.
    this->resourceSpineData.push_back(std::make_unique<ResourceSynapseSpine>());
    ResourceSynapseSpine* newSpine = resourceSpineData.back().get();
    //Step weights has been removed fron here
    //newSynapse->SetWeight(this->GenerateSynapticWeight());

    //this->weightsSum += newSynapse->GetWeight();
    newSpine->SetIdInMorpho(this->spineIdGenerator++);
    //Branch
    int branch {AllocateBranch(branchTarget)};
    int position{PopSynapseSlotFromBranch(branch, branchTarget.firstSlotTrueLastSlotFalse)};
    newSpine->SetBranchPositionId(position);
    newSpine->SetBranchId(branch);
    newSpine->SetDistanceFromNode(position*branches.at(branch)->synapticGap);//This has to be updated if we switch to double 
    newSpine->SetAlphaStimBump(baseAlphaStimBump);
    newSpine->SetBiasLTD(biasLTD);
    newSpine->SetAlphaBasal(alphaBasal);
    newSpine->SetAlphaExpDecay(alphaStimulusExpDecay);
    branches.at(branch)->synapseSlotClosedIndex.push_back(position);
    //branches.at(branch)->morphoSynapseIDs.push_back(newSynapse->GetIdInMorpho());
    //branches.at(branch)->synapseSlotToMorphoIDMap[position]=newSynapse->GetIdInMorpho();
    //Storage (other)
    this->baseSpineData.push_back(static_cast<BaseSpinePtr>(newSpine));
    this->branchedSpineData.push_back(static_cast<BranchedSpinePtr>(newSpine));


    return this->baseSpineData.back();
}

std::vector<double> TraceRBranchedHSTDP::GetOverallSynapticProfile() const {
    /*
     * returned array organised as follows:
     * item 1: average synaptic weight
     * item 2: totalLTD Events
     * item 3: totalLTP Events
     * item 4: average plasticity events
     * */
    std::vector<double> dataArray(3);
    size_t noSpines{this->baseSpineData.size()};
    //size_t sizeOfSpineData {this->baseSpineData.size()};
    // double weightSum = std::accumulate(this->baseSpineData.begin(), this->baseSpineData.end(), 0.0, [] (double accumulator, const BaseSpinePtr syn) {
    //                                     return accumulator + syn->GetWeight(); 
    //                                     });
    //CalcMorphoPlasticityEvents();

   dataArray.at(0) = std::accumulate(this->baseSpineData.begin(), this->baseSpineData.end(), 0.0, [] (double accumulator, const BaseSpinePtr syn) {
       return accumulator + syn->GetWeightUncoupled();}) / noSpines;
   dataArray.at(1) = this->totalPostSpikes;
   dataArray.at(2) = this->totalPreSpikes;
   return dataArray;
}

std::string TraceRBranchedHSTDP::GetOverallSynapticProfileHeaderInfo() const {
    return std::string("{<average weight>, <total post spikes> ,<total pre spikes>}");
}

// void BranchedResourceHeteroSTDP::CalcMorphoPlasticityEvents()
// {
//     totalLTDEvents = std::accumulate(this->branches.begin(), this->branches.end(), 0,
//                                        [] (int accumulator, const BranchPtr branch) { return accumulator + branch->LTDevents; });
//     totalLTPEvents = std::accumulate(this->branches.begin(), this->branches.end(), 0,
//                                        [] (int accumulator, const BranchPtr branch) { return accumulator + branch->LTPevents; });
//     totalPlasticityEvents = totalLTDEvents + totalLTPEvents;
// }
