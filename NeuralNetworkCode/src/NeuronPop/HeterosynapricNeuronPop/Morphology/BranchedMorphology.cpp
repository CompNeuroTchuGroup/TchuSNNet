#include "BranchedMorphology.hpp"



BranchedMorphology::BranchedMorphology(GlobalSimInfo * info): Morphology(info){

}

void BranchedMorphology::recordPostSpike() {
    Morphology::recordPostSpike();
    std::fill(this->integratePostSpike.begin(), this->integratePostSpike.end(), true);
    this->postSpiked = true;
}

void BranchedMorphology::recordExcitatoryPreSpike(unsigned long synSpikerId) {
    Morphology::recordExcitatoryPreSpike(synSpikerId);
    this->integratePreSpike.at(synSpikerId) = true;
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
            this->synapticGap = std::stoi(values.at(0));
            synapticGapInitialized = true;
        } else if (name.find("distribute_weights") != std::string::npos) {
            //This whole part is experimental, it seems it was not completely tested
            //As such, this is deprecated from publication
            if (values.at(0) == "true") {
                distributeWeights = true;
            } else {
                    this->initialWeights = std::stod(values.at(1));
            }
        } else if (name.find("dendrite_branchings") != std::string::npos) {
            this->branchings = std::stoi(values.at(0));
            if (this->branchings>28){
                //EXCEPTION, integer overflow.
            }
            branchingsInitialized=true;
        }
    }

    assertm(dendriteInitialized, "Using heterosynaptic synapses without specifying dendritic_length is not allowed.");
    assertm(synapticGapInitialized, "Using heterosynaptic synapses without specifying synaptic_gap is not allowed.");
    assertm(branchingsInitialized, "Using branched morphology with no branchings specified.");
    setUpBranchings(this->branchings);
    for (auto branch : this->branches){
        setUpSynapseSlots(branch);
    }
}

void BranchedMorphology::SaveParameters(std::ofstream *stream, std::string neuronPreId) {
    Morphology::SaveParameters(stream, neuronPreId);

    *stream << neuronPreId<<"_morphology_branch_length\t\t"<<std::to_string(this->branchLength)<<" #μm";
    *stream << "\t"<<"#Length of each branch.\n";

    *stream << neuronPreId<<"_morphology_synaptic_gap\t\t\t"<<std::to_string(this->synapticGap)<<" #μm";
    *stream << "\t"<<"#Unitary distance between synapse slots.\n";

    *stream << neuronPreId<<"_morphology_distribute_weights\t\t";
    if (this->distributeWeights){
        *stream << "true\t";
    }else{
        *stream<<"false\t"<<std::to_string(this->initialWeights);
    }
    *stream << "\t"<<"#The bool corresponds to distributing weight between min and max uniformally. The number will be the weight assigned to all synapses if bool is false (do not confuse with implementation in MonoDendriteSTDP).\n";
    *stream << neuronPreId<<"_morphology_dendrite_branchings\t\t"<<std::to_string(this->branchings);
    *stream << "\t"<<"#This specifies the number of branchings in the dendritic tree FOR EVERY EXISTING BRANCH. Total isolated branches are 2^n. More than 28 will cause integer overflow\n";
}


double BranchedMorphology::generateSynapticWeight(){
    double weight{};
    std::uniform_real_distribution<double> distribution(this->minWeight,this->maxWeight);
            if (distributeWeights) {
                std::default_random_engine& generator = this->info->globalGenerator;
                weight = distribution(generator);
            } else {
                weight = this->initialWeights; // assuming a range of weight between 0 and 2, weight is initialized to midpoint: 1
            }
        this->weightsSum += weight;
        return weight;
}

int BranchedMorphology::randomBranchAllocation()
{
        std::default_random_engine& generator = this->info->globalGenerator;
        //For now, the distribution will be uniform
        std::uniform_int_distribution<int> branchdsitribution(0,static_cast<int>(branches.size())+1);
        int branchID{branchdsitribution(generator)};
        return branchID;
}

int BranchedMorphology::orderedBranchAllocation()
{
    for (auto& branch:branches){
        if (branch->openSynapsesSlots.size()==0){
            continue;
        }
        return branch->branchId;
    }
}

void BranchedMorphology::RandomSynapseAllocation(std::shared_ptr<Branch> branch)
{
    std::default_random_engine& generator = this->info->globalGenerator;
        for (auto& branch:branches){
            std::vector<int> possibleSlots(branch->spikedSyn.size());
            std::iota(possibleSlots.begin(), possibleSlots.end(), 0);
            //Now we have our vector from 0 to maxSlots to pull random numbers from
            std::sample(possibleSlots.begin(), possibleSlots.end(),std::back_inserter(branch->openSynapsesSlots),branch->spikedSyn.size(),this->info->globalGenerator);
        //Then I will have to pop_front() in allocateNewSynapse
        }
}

void BranchedMorphology::OrderedSynapseAllocation(std::shared_ptr<Branch> branch)
{
    for (auto& branch:branches){
        std::deque<int> possibleSlots(branch->spikedSyn.size());
        std::iota(possibleSlots.begin(), possibleSlots.end(), 0);
        //Now we have our vector from 0 to maxSlots to pull random numbers from
        copy(possibleSlots.begin(), possibleSlots.end(), back_inserter(branch->openSynapsesSlots));
    }
    return;
}

/*void BranchedMorphology::AlternatedSynapseAllocation(std::shared_ptr<Branch> branch)
{
        for (auto& branch:branches){
        std::deque<int> possibleSlots(branch->spikedSyn.size());
        std::iota(possibleSlots.begin(), possibleSlots.end(), 0);
        copy(possibleSlots.begin(), possibleSlots.end(), back_inserter(branch->openSynapsesSlots));
        }
        //For now it is the same as ordered, the possibility of alternating synapses will probably be more useful with multiple synapses pre-post neuron.
}*/

std::valarray<double> BranchedMorphology::getIndividualSynapticProfile(unsigned long synapseId) const {
    /*
     * returned array organised as follows:
     * item 1: distance of synapse from branch root
     * item 2: value of the synaptic weight
     * item 3: last spike time of the synapse
     * */
    return synapseData.at(synapseId)->getIndividualSynapticProfile();
}

std::valarray<double> BranchedMorphology::getOverallSynapticProfile() const {
    /*
     * returned array organised as follows:
     * item 1: average synaptic weight
     * item 2: total post spikes
     * item 3: total pre spikes
     * */
    std::valarray<double> ret(3);

    double weightSum = std::accumulate(this->synapseData.begin(), this->synapseData.end(), 0.0,
                                       [] (const double acc, const std::shared_ptr<SynapseSpine>& syn) { return acc + syn->getWeight(); });

    ret[0] = weightSum / this->synapseData.size();
    ret[1] = this->totalPostSpikes;
    ret[2] = this->totalPreSpikes;
    return ret;
}
void BranchedMorphology::setUpBranchings (int remainingBranchingEvents, std::vector<int> anteriorBranches){ 
    //This is a recursive function that sets up the branched dendritic tree and is generalized for 0 branchings (1 branch). This function has been unit tested by Antoni.
    remainingBranchingEvents-=1;
    //First call is done with an empty int vector
    for (int i : {1,2}){
        int branchId{this->generateBranchId()};
        this->branches.emplace_back(std::make_shared<Branch>(Branch(this->synapticGap, this->branchLength, anteriorBranches, branchId)));//This vector should be sorted by ID by default (tested).
        //Constructor here
        if(remainingBranchingEvents>0){
            std::vector<int> anteriorBranchesCopy(anteriorBranches);
            anteriorBranchesCopy.push_back(branchId);
            this->setUpBranchings(remainingBranchingEvents, anteriorBranchesCopy);
        }
    }
}