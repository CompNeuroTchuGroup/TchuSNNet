//
// Created by Antoni Bertolin on 14.06.23
//
#include "./BranchedMorphology.hpp"
#include "BranchedMorphology.hpp"
BranchedMorphology::BranchedMorphology(GlobalSimInfo* infoGlobal): Morphology(infoGlobal){

}

void BranchedMorphology::RecordPostSpike() {
    Morphology::RecordPostSpike();
}

void BranchedMorphology::RecordExcitatoryPreSpike(int spikedSynapseId) {
    //This function is NOT DELAY COMPATIBLE (careful with the delays in synapse objects)
    Morphology::RecordExcitatoryPreSpike(spikedSynapseId);
    // this->branches.at(this->branchedSpineData.at(spikedSynapseId)->GetBranchId())->spikedSyn.at(this->branchedSpineData.at(spikedSynapseId)->GetBranchPositionId())=true;
    BranchedSpinePtr synapseSpine = this->branchedSpineData.at(spikedSynapseId);
    this->branches.at(synapseSpine->GetBranchId())->spikedSpinesInTheBranch.push_back(synapseSpine->GetBranchPositionId());
}

void BranchedMorphology::Reset() {
    std::for_each(branches.begin(), branches.end(), [] (BranchPtr branch){
        branch->spikedSpinesInTheBranch.clear();
    });
    this->postSpiked=false;
}

void BranchedMorphology::LoadParameters(const std::vector<FileEntry>& morphologyParameters) {
    Morphology::LoadParameters(morphologyParameters);

    bool branchLengthInitiallized = false,
        synapticGapInitialized = false,
        branchingsInitialized = false;

    for (auto& [parameterName, parameterValues] : morphologyParameters) {
        if(parameterName.find("branchLength") != std::string::npos){
            this->branchLength = std::stod(parameterValues.at(0));
            branchLengthInitiallized = true;
        } else if (parameterName.find("synapticGap") != std::string::npos) {
            this->synapticGap = std::stod(parameterValues.at(0));
            synapticGapInitialized = true;
        } else if (parameterName.find("distributeWeights") != std::string::npos) {
            //This whole part is experimental, it seems it was not completely tested
            //As such, this is deprecated from publication
            if (parameterValues.at(0) == "true") {
                this->distributeWeights = true;
            } else if (parameterValues.at(0) == "false") {
                this->distributeWeights = false;
                    // this->initialWeights = std::stod(parameterValues.at(1));
            } else {
                throw "Logical error in distributeWeights";
            }
        } else if (parameterName.find("dendriteBranchings") != std::string::npos) {
            this->branchings = std::stoi(parameterValues.at(0));
            if (this->branchings>28){
                //EXCEPTION, integer overflow.
                throw "Dendrite branching causing integer overflow";
            }
            branchingsInitialized=true;
        } else if (parameterName.find("synapseAllocation") != std::string::npos){
            if (parameterValues.at(0).find("ordered") != std::string::npos){
                this->orderedSpineAllocationB=true;
            }else if (parameterValues.at(0).find("random") != std::string::npos){
                this->randomSpineAllocationB=true;
            }
        } /*else if (parameterName.find("branch_allocation") != std::string::npos){
            if (parameterValues.at(0).find("ordered") != std::string::npos){
                OrderedBranchAllocationB=true;
            } else if (parameterValues.at(0).find("set") != std::string::npos){
                setBranchAllocationB=true;
            } else if (parameterValues.at(0).find("random") != std::string::npos){
                RandomBranchAllocationB=true;
            }
        }*/
    }
    //Improve exception handling
    if (!branchLengthInitiallized || !synapticGapInitialized||!branchingsInitialized){throw "Variable not initialized";}
    if (!orderedSpineAllocationB && !randomSpineAllocationB){throw "Logical error in BranchedMorphology";}
    //if (!OrderedBranchAllocationB && !RandomBranchAllocationB && !setBranchAllocationB){throw;}
    SetUpBranchedMorphology();
}

void BranchedMorphology::CheckParameters(const std::vector<FileEntry> &parameters) {
    Morphology::CheckParameters(parameters);
    for (auto& [parameterName, parameterValues] : parameters) {
        if(parameterName.find("branchLength") != std::string::npos){
            if(!(this->branchLength == std::stod(parameterValues.at(0)))){
                throw "Branch length was not consistent in plasticity model parameters.";
            }
        } else if (parameterName.find("synapticGap") != std::string::npos) {
            if(!(this->synapticGap == std::stod(parameterValues.at(0)))){
                throw "Branch length was not consistent in plasticity model parameters.";
            }
        } else if (parameterName.find("distributeWeights") != std::string::npos) {
            //This whole part is experimental, it seems it was not completely tested
            //As such, this is deprecated from publication
            if (parameterValues.at(0) == "true") {
                if(!(this->distributeWeights)){
                    throw"distributeWeights was not consistent in plasticity model parameters.";
                }
            } else if (parameterValues.at(0) == "false") {
                if(this->distributeWeights){
                    throw "distributeWeights was not consistent in plasticity model parameters.";
                }
            }
        } else if (parameterName.find("dendriteBranchings") != std::string::npos) {
            if(this->branchings != std::stoi(parameterValues.at(0))){
                throw "branchings was not consistent in plasticity model parameters.";
            }
        } else if (parameterName.find("synapseAllocation") != std::string::npos){
            if (parameterValues.at(0).find("ordered") != std::string::npos){
                if(!(this->orderedSpineAllocationB)){
                    throw "synapseAllocation was not consistent in plasticity model parameters.";      
                }
            } else if (parameterValues.at(0).find("random") != std::string::npos){
                if(!(this->randomSpineAllocationB)){
                    throw "synapseAllocation was not consistent in plasticity model parameters.";      
                }
            }
        }
    }
}
void BranchedMorphology::SetUpBranchedMorphology() {
    //If a different starting point for binary branch splits is desired, add more branches here and remember to add them in the recursive function call.
    SetUpBranchings(this->branchings);
    std::for_each(branches.begin(), branches.end(), [this] (BranchPtr const branch){
        SetUpSynapseSlots(branch);
    });
}

void BranchedMorphology::SetUpSynapseSlots(BranchPtr branch) {
    if (this->orderedSpineAllocationB){
        this->OrderedSynapseAllocation(branch);
    } else if (this->randomSpineAllocationB){
        this->RandomSynapseAllocation(branch);
    } else {
        throw "Logic error in BranchedMorphology";
    }
}

void BranchedMorphology::SaveParameters(std::ofstream& wParameterFile, std::string neuronIdentificator) const {
    Morphology::SaveParameters(wParameterFile, neuronIdentificator);

    wParameterFile << neuronIdentificator<<"branchLength\t\t"<<std::to_string(this->branchLength)<<" #μm";
    wParameterFile << "\t"<<"#Length of each branch.\n";

    wParameterFile << neuronIdentificator<<"synapticGap\t\t\t"<<std::to_string(this->synapticGap)<<" #μm";
    wParameterFile << "\t"<<"#Distance between synapse spines.\n";

    wParameterFile << neuronIdentificator<<"distributeWeights\t\t"<<std::boolalpha<<this->distributeWeights<<std::noboolalpha<<" #μm";
    wParameterFile << "\t"<<"#Distance between synapse spines.\n";

    // wParameterStream << neuronIdentificator<<"distributeWeights\t\t";
    // if (this->distributeWeights){
    //     wParameterStream << "true\t";
    // }else{
    //     wParameterStream<<"false\t"<<std::to_string(this->initialWeights);
    // }
    // wParameterStream << "\t"<<"#The bool corresponds to distributing weight between min and max uniformally. The number will be the weight assigned to all synapses if bool is false (do not confuse with implementation in MonoDendriteSTDP).\n";
    
    wParameterFile << neuronIdentificator<<"dendriteBranchings\t\t"<<std::to_string(this->branchings);
    wParameterFile << "\t\t"<<"#This specifies the number of branchings in the dendritic tree FOR EVERY EXISTING BRANCH. Total isolated branches are 2^n. More than 28 will cause integer overflow\n";
    
    wParameterFile << neuronIdentificator<<"synapseAllocation\t\t";
    if (this->orderedSpineAllocationB){
        wParameterFile << "ordered\t";
    }else if (this->randomSpineAllocationB){
        wParameterFile<<"random\t";
    }
    wParameterFile << "\t"<<"#'ordered' synapse allocation will allocate synapses from the branch node to the end of the branch. 'random' will allocate random positions in each branch\n";

    //wParameterStream << neuronIdentificator<<"seed\t\t\t"<<std::to_string(this->seed)<<"\n";//Missing comments
    /*
    wParameterStream << neuronIdentificator<<"branch_allocation\t\t";
    if (this->RandomBranchAllocationB){
        wParameterStream << "random\t";
    }else if (this->OrderedBranchAllocationB){
        wParameterStream<<"ordered\t";
    }else if (this->setBranchAllocationB){
        wParameterStream<<"set\t";
    }
    wParameterStream << "\t"<<"#Missing comments\n";*/
}


// BaseSpinePtr BranchedMorphology::AllocateNewSynapse(const BranchTargeting& branchTarget){
//     BranchedSpinePtr newSynapse = std::make_shared<BranchedSynapseSpine>();

//     //REFORMAT, REWRITE WITH CONSTRUCTOR    
//     //Step weights has been removed fron here
//     newSynapse->SetWeight(this->GenerateSynapticWeight());

//     //this->weightsSum += newSynapse->GetWeight();
//     newSynapse->SetIdInMorpho(this->spineIdGenerator++);
//     //Branch
//     int branch {AllocateBranch(branchTarget)};
//     if (branches.at(branch)->openSpineSlots.empty()){
//         throw "No synapses available on dendrite for new allocation.";
//     }
//     //Position
//     int position{branches.at(branch)->openSpineSlots.front()};
//     branches.at(branch)->openSpineSlots.pop_front();
//     newSynapse->SetBranchPositionId(position);
//     newSynapse->SetDistanceFromNode(position*branches.at(branch)->synapticGap);
//     branches.at(branch)->synapseSlotClosedIndex.push_back(position);
//     //branches.at(branch)->morphoSynapseIDs.push_back(newSynapse->GetIdInMorpho());
//     //branches.at(branch)->synapseSlotToMorphoIDMap.at(position)=newSynapse->GetIdInMorpho();
//     //Storage (other)
//     this->baseSpineData.push_back(static_pointer_cast<BaseSynapseSpine>(newSynapse));
//     this->branchedSpineData.push_back(newSynapse);

//     //this->spikedSynapses.push_back(false);
//     //this->integratePostSpike.push_back(false);
//     //this->integratePreSpike.push_back(false);

//     return static_pointer_cast<BaseSynapseSpine>(newSynapse);
// }

int BranchedMorphology::AllocateBranch(const BranchTargeting branchTarget) {
        if (branchTarget.setTargetBranch && !(branchTarget.targetBranch>=static_cast<int>(this->branches.size()))){//SHould this condition be equal or less than?
            return branchTarget.targetBranch;
        } else if (branchTarget.randomTargetBranch){
            return RandomBranchAllocation();
        } else if (branchTarget.orderedTargetBranch){
            return OrderedBranchAllocation();
        } else {
            throw "There was no synaptic targeting detected.";
        }
}

int BranchedMorphology::RandomBranchAllocation() {
        //For now, the distribution will be uniform
        std::uniform_int_distribution<int> branchDistribution(0,static_cast<int>(branches.size()-1));
        int branchID{branchDistribution(this->generator)};
        if (this->branches.at(branchID)->openSpineSlots.empty()){
            branchID=RandomBranchAllocation();
        }
        return branchID;
}

int BranchedMorphology::OrderedBranchAllocation() {
    int branchId{};
    for (BranchPtr branch : branches){
        if (!branch->openSpineSlots.empty()){
            branchId = branch->branchId ;
            break;
        }
    }
    return branchId;
}

void BranchedMorphology::RandomSynapseAllocation(BranchPtr branch) {
    //std::mt19937& generator = this->generator;
    std::vector<int> possibleSlots(branch->branchSlots);
    std::iota(possibleSlots.begin(), possibleSlots.end(), 0);
   //Now we have our vector from 0 to maxSlots to pull random numbers from
    std::sample(possibleSlots.begin(), possibleSlots.end(),std::back_inserter(branch->openSpineSlots),branch->branchSlots,generator);
    //Then I will have to pop_front() in AllocateNewSynapse
}

void BranchedMorphology::OrderedSynapseAllocation(BranchPtr branch) {
    std::deque<int> possibleSlots(branch->branchSlots);
    std::iota(possibleSlots.begin(), possibleSlots.end(), 0);
    //Now we have our vector from 0 to maxSlots to pull random numbers from
    copy(possibleSlots.begin(), possibleSlots.end(), back_inserter(branch->openSpineSlots));
    //Then I will have to pop_front() in AllocateNewSynapse
}


void BranchedMorphology::SetUpBranchings(int remainingBranchingEvents, std::vector<int> anteriorBranches) {
    //This is a recursive function that sets up the branched dendritic tree and is generalized for 0 branchings (1 branch). This function has been unit tested by Antoni.
    remainingBranchingEvents-=1;
    //First call is done with an empty int vector
    for (int iteration = 0; iteration < 2;iteration++) {
        int branchId{this->GenerateBranchId()};
        this->branches.push_back(new Branch(this->synapticGap, this->branchLength, anteriorBranches, branchId));//This vector should be sorted by ID by default (tested).
        //Constructor here
        if(remainingBranchingEvents>0){
            std::vector<int> anteriorBranchesCopy(anteriorBranches);
            anteriorBranchesCopy.push_back(branchId);
            this->SetUpBranchings(remainingBranchingEvents, anteriorBranchesCopy);
        }
    }
    //Here we can create the node pointers:
}

void BranchedMorphology::PostConnectSetUp()
{
    //Here we do all the function calls that could not be done in the constructor/LP. 
    //This is done to adapt to the fact that synapses do not exist until ConnectNeurons() is called in the NeuralNetwork::Simulate() function.
    std::for_each(branches.begin(), branches.end(), [this] (BranchPtr branch){
        branch->PostConnectSetUp(branchedSpineData);
    });
}
