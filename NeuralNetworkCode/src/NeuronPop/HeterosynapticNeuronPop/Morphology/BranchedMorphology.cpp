#include "./BranchedMorphology.hpp"
#include "BranchedMorphology.hpp"

BranchedMorphology::BranchedMorphology(GlobalSimInfo * info): Morphology(info){

}

void BranchedMorphology::RecordPostSpike() {
    Morphology::RecordPostSpike();
}

void BranchedMorphology::RecordExcitatoryPreSpike(int spikedSynapseId) {
    //This function is NOT DELAY COMPATIBLE (careful with the delays in synapse objects)
    Morphology::RecordExcitatoryPreSpike(spikedSynapseId);
    // this->branches.at(this->branchedSynapseData.at(spikedSynapseId)->GetBranchId())->spikedSyn.at(this->branchedSynapseData.at(spikedSynapseId)->GetBranchPositionId())=true;
    BranchedSpinePtr& synapseSpine = this->branchedSynapseData.at(spikedSynapseId);
    this->branches.at(synapseSpine->GetBranchId())->spikedSynapsesInTheBranch.push_back(synapseSpine->GetBranchPositionId());
}

void BranchedMorphology::Reset()
{
    for (BranchPtr branch : branches){
        // std::fill(branch->spikedSyn.begin(), branch->spikedSyn.end(), false);
        branch->spikedSynapsesInTheBranch.clear();
    }
    this->postSpiked=false;
}

void BranchedMorphology::LoadParameters(std::vector<std::string> *input) {
    Morphology::LoadParameters(input);

    std::string name;
    std::vector<std::string> values;

    bool dendriteInitialized = false,
        synapticGapInitialized = false,
        branchingsInitialized = false;

    for (auto & it : *input) {
        SplitString(&it,&name,&values);

        if(name.find("branch_length") != std::string::npos){
            this->branchLength = std::stoi(values.at(0));
            dendriteInitialized = true;
        } else if (name.find("synaptic_gap") != std::string::npos) {
            this->synapticGap = std::stod(values.at(0));
            synapticGapInitialized = true;
        } else if (name.find("distribute_weights") != std::string::npos) {
            //This whole part is experimental, it seems it was not completely tested
            //As such, this is deprecated from publication
            if (values.at(0) == "true") {
                this->distributeWeights = true;
            } else {
                    this->initialWeights = std::stod(values.at(1));
            }
        } else if (name.find("dendrite_branchings") != std::string::npos) {
            this->branchings = std::stoi(values.at(0));
            if (this->branchings>28){
                //EXCEPTION, integer overflow.
                throw;
            }
            branchingsInitialized=true;
        } else if (name.find("synapse_allocation") != std::string::npos){
            if (values.at(0).find("ordered") != std::string::npos){
                this->orderedSynAllocationB=true;
            }else if (values.at(0).find("random") != std::string::npos){
                this->randomSynAllocationB=true;
            }
        } /*else if (name.find("branch_allocation") != std::string::npos){
            if (values.at(0).find("ordered") != std::string::npos){
                OrderedBranchAllocationB=true;
            } else if (values.at(0).find("set") != std::string::npos){
                setBranchAllocationB=true;
            } else if (values.at(0).find("random") != std::string::npos){
                RandomBranchAllocationB=true;
            }
        }*/
    }
    //Improve exception handling
    assertm(dendriteInitialized, "Using heterosynaptic synapses without specifying dendritic_length is not allowed.");
    assertm(synapticGapInitialized, "Using heterosynaptic synapses without specifying synaptic_gap is not allowed.");
    assertm(branchingsInitialized, "Using branched morphology with no branchings specified.");
    if (!orderedSynAllocationB && !randomSynAllocationB){throw;}
    //if (!OrderedBranchAllocationB && !RandomBranchAllocationB && !setBranchAllocationB){throw;}
    SetUpBranchedMorphology();
}

void BranchedMorphology::SetUpBranchedMorphology()
{
    SetUpBranchings(this->branchings);
    for (auto& branch : this->branches){
        SetUpSynapseSlots(branch);
    }
}

void BranchedMorphology::SetUpSynapseSlots(BranchPtr branch)
{
    if (this->orderedSynAllocationB){
        this->OrderedSynapseAllocation(branch);
    } else if (this->randomSynAllocationB){
        this->RandomSynapseAllocation(branch);
    }
}

void BranchedMorphology::SaveParameters(std::ofstream *stream, std::string neuronPreId) {
    Morphology::SaveParameters(stream, neuronPreId);

    *stream << neuronPreId<<"_morphology_branch_length\t\t"<<std::to_string(this->branchLength)<<" #μm";
    *stream << "\t"<<"#Length of each branch.\n";

    *stream << neuronPreId<<"_morphology_synaptic_gap\t\t\t"<<std::to_string(this->synapticGap)<<" #μm";
    *stream << "\t"<<"#Distance between synapse spines.\n";

    // *stream << neuronPreId<<"_morphology_distribute_weights\t\t";
    // if (this->distributeWeights){
    //     *stream << "true\t";
    // }else{
    //     *stream<<"false\t"<<std::to_string(this->initialWeights);
    // }
    // *stream << "\t"<<"#The bool corresponds to distributing weight between min and max uniformally. The number will be the weight assigned to all synapses if bool is false (do not confuse with implementation in MonoDendriteSTDP).\n";
    
    *stream << neuronPreId<<"_morphology_dendrite_branchings\t\t"<<std::to_string(this->branchings);
    *stream << "\t"<<"#This specifies the number of branchings in the dendritic tree FOR EVERY EXISTING BRANCH. Total isolated branches are 2^n. More than 28 will cause integer overflow\n";
    
    *stream << neuronPreId<<"_morphology_synapse_allocation\t\t";
    if (this->orderedSynAllocationB){
        *stream << "ordered\t";
    }else if (this->randomSynAllocationB){
        *stream<<"random\t";
    }
    *stream << "\t"<<"#'ordered' synapse allocation will allocate synapses from the branch node to the end of the branch. 'random' will allocate random positions in each branch\n";

    //*stream << neuronPreId<<"_morphology_seed\t\t\t"<<std::to_string(this->seed)<<"\n";//Missing comments
    /*
    *stream << neuronPreId<<"_morphology_branch_allocation\t\t";
    if (this->RandomBranchAllocationB){
        *stream << "random\t";
    }else if (this->OrderedBranchAllocationB){
        *stream<<"ordered\t";
    }else if (this->setBranchAllocationB){
        *stream<<"set\t";
    }
    *stream << "\t"<<"#Missing comments\n";*/
}


BaseSpinePtr BranchedMorphology::AllocateNewSynapse(const HeteroCurrentSynapse& synapse)//This should not work
{
    BranchedSpinePtr newSynapse = std::make_shared<BranchedSynapseSpine>();

    //REFORMAT, REWRITE WITH CONSTRUCTOR    
    //Step weights has been removed fron here
    newSynapse->SetWeight(this->GenerateSynapticWeight());

    //this->weightsSum += newSynapse->GetWeight();
    newSynapse->SetIdInMorpho(this->synapseIdGenerator++);
    //Branch
    int branch {AllocateBranch(synapse)};
    if (branches.at(branch)->openSynapsesSlots.size()==0){
        throw noAllocatableSynapseException();
    }
    //Position
    int position{branches.at(branch)->openSynapsesSlots.front()};
    branches.at(branch)->openSynapsesSlots.pop_front();
    newSynapse->SetBranchPositionId(position);
    newSynapse->SetDistanceFromNode(position*branches.at(branch)->synapticGap);
    branches.at(branch)->synapseSlotClosedIndex.push_back(position);
    //branches.at(branch)->morphoSynapseIDs.push_back(newSynapse->GetIdInMorpho());
    branches.at(branch)->synapseSlotToMorphoIDMap.at(position)=newSynapse->GetIdInMorpho();
    //Storage (other)
    this->baseSynapseData.push_back(static_cast<BaseSpinePtr>(newSynapse));
    this->branchedSynapseData.push_back(newSynapse);

    //this->spikedSynapses.push_back(false);
    //this->integratePostSpike.push_back(false);
    //this->integratePreSpike.push_back(false);

    return static_cast<BaseSpinePtr>(newSynapse);
}

int BranchedMorphology::AllocateBranch(const HeteroCurrentSynapse &synapse)
{
        const BranchTargeting& target= synapse.getBranchTarget();
        if (target.setTargetBranch && !(target.targetBranch>=this->branches.size())){//SHould this condition be equal or less than?
            return target.targetBranch;
        } else if (target.randomTargetBranch){
            return RandomBranchAllocation();
        } else if (target.orderedTargetBranch){
            return OrderedBranchAllocation();
        } else {
            throw noAllocatableSynapseException();
        }
}

int BranchedMorphology::RandomBranchAllocation()
{
        std::mt19937& generator = this->generator;
        //For now, the distribution will be uniform
        std::uniform_int_distribution<int> branchDistribution(0,static_cast<int>(branches.size()-1));
        int branchID{branchDistribution(generator)};
        if (this->branches.at(branchID)->openSynapsesSlots.size()==0){
            branchID=RandomBranchAllocation();
        }
        return branchID;
}

int BranchedMorphology::OrderedBranchAllocation()
{
    int branchId{};
    for (auto& branch:branches){
        if (branch->openSynapsesSlots.size()==0){
            continue;
        }
        branchId= branch->branchId ;
        break;
    }
    return branchId;
}

void BranchedMorphology::RandomSynapseAllocation(BranchPtr& branch)
{
    //std::mt19937& generator = this->generator;
    std::vector<int> possibleSlots(branch->branchSlots);
    std::iota(possibleSlots.begin(), possibleSlots.end(), 0);
   //Now we have our vector from 0 to maxSlots to pull random numbers from
    std::sample(possibleSlots.begin(), possibleSlots.end(),std::back_inserter(branch->openSynapsesSlots),branch->branchSlots,this->info->globalGenerator);
    //Then I will have to pop_front() in AllocateNewSynapse
}

void BranchedMorphology::OrderedSynapseAllocation(BranchPtr& branch)
{
    std::deque<int> possibleSlots(branch->branchSlots);
    std::iota(possibleSlots.begin(), possibleSlots.end(), 0);
    //Now we have our vector from 0 to maxSlots to pull random numbers from
    copy(possibleSlots.begin(), possibleSlots.end(), back_inserter(branch->openSynapsesSlots));
    //Then I will have to pop_front() in AllocateNewSynapse
}

/*void BranchedMorphology::AlternatedSynapseAllocation(BranchPtr branch)
{
        for (auto& branch:branches){
        std::deque<int> possibleSlots(branch->branchSlots);
        std::iota(possibleSlots.begin(), possibleSlots.end(), 0);
        copy(possibleSlots.begin(), possibleSlots.end(), back_inserter(branch->openSynapsesSlots));
        }
        //For now it is the same as ordered, the possibility of alternating synapses will probably be more useful with multiple synapses pre-post neuron.
}*/
// void BranchedMorphology::CalcMorphoPlasticityEvents()
// {
//     totalPlasticityEvents=std::accumulate(this->branches.begin(), this->branches.end(), 0,//UNRESOLVED, does this give intended output?
//                                        [] (unsigned long acc, const BranchPtr& branch) { return acc + branch->LTPevents + branch->LTDevents; });
// }

void BranchedMorphology::SetUpBranchings(int remainingBranchingEvents, std::vector<int> anteriorBranches)
{
    //This is a recursive function that sets up the branched dendritic tree and is generalized for 0 branchings (1 branch). This function has been unit tested by Antoni.
    remainingBranchingEvents-=1;
    //First call is done with an empty int vector
    for (int i = 0; i < 2;i++) {
        int branchId{this->GenerateBranchId()};
        this->branches.emplace_back(std::make_shared<Branch>(Branch(this->synapticGap, this->branchLength, anteriorBranches, branchId)));//This vector should be sorted by ID by default (tested).
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
    for (BranchPtr branch: branches){
        branch->PostConnectSetUp(branchedSynapseData);
    }
}
