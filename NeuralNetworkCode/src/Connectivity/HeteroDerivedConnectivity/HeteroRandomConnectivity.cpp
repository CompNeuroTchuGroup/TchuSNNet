#include "./HeteroRandomConnectivity.hpp"
#include "../../Synapse/Synapse.hpp"
#include "../../Synapse/HeteroCurrentSynapse.hpp"

HeteroRandomConnectivity::HeteroRandomConnectivity(Synapse* syn, GlobalSimInfo* info):
    RandomConnectivity(syn, info), HeteroConnectivity(syn->GetNoNeuronsPre()) {
}

void HeteroRandomConnectivity::ConnectNeurons() {

    unsigned long numPostNeurons = this->synapse->GetNoNeuronsPost();

    unsigned long output_Interval = numPostNeurons / 10;
    output_Interval = output_Interval == 0 ? 1 : output_Interval; //Comparing a long to ints is not going to work out. Must be changed.

    std::uniform_int_distribution<int> distribution(0,(int)this->synapse->GetNoNeuronsPre()-1);

    unsigned long source, countedSourceNeurons, synapseId;

    for(unsigned long target = 0; target < numPostNeurons; target++) {
        countedSourceNeurons = 0;

        while(countedSourceNeurons < this->noSourceNeurons) {// This works in recurrent synapses because the noSourceNeurons is -1
            source = distribution(generator);

            if(this->synapse->IsRecurrent() && source == target) {
                continue;
            }

            if(!target_id[source].empty() && target_id[source].back() == target) {//Has it been already plugged? (because of ordering)
                continue;
            }

            try {
                synapseId = this->synapse->allocateSynapse(source, target);// synID equal to id in HCS
                synapticTargets[source].push_back(std::make_pair(target, synapseId)); //Multi targeting has to change here (for random connectivity)
                target_id[source].push_back(target);//Are both truly necessary?
                countedSourceNeurons++;
            } catch (...) {//const noAllocatableSynapseException& e
                throw;
            }
        }
        if ((target % output_Interval) == 0) {
            std::cout << 100 * target / numPostNeurons << "%-";
        }
    }
    if (this->synapse->GetNeuronsPost()->IsBranchedBool()){
        this->synapse->GetNeuronsPost()->PostConnectSetUp();
    }
    std::cout << "100%\n";

}

const std::string HeteroRandomConnectivity::GetTypeStr() {
    return str_heteroRandomConnectivity;
}

std::vector<std::pair<unsigned long, unsigned long>>& HeteroRandomConnectivity::getSynapticTargets(const unsigned long& preId){
    return HeteroConnectivity::getSynapticTargets(preId);
}


//unsigned long HeteroRandomConnectivity::GetNumberAverageSourceNeurons() {
//    return 0;
//}
//
//int* HeteroRandomConnectivity::GetDistributionD(long preNeuronId, long postNeuronId) {
//    return &D_distribution[preNeuronId][postNeuronId];
//}
//
//void HeteroRandomConnectivity::SetDistributionD() {
//    this->D_distribution->resize(this->synapse->GetNoNeuronsPre());
//    for (unsigned long i = 0; i < this->synapse->GetNoNeuronsPre(); ++i) {
//        this->D_distribution[i].resize(this->synapse->GetNoNeuronsPost());
//        std::fill(this->D_distribution[i].begin(), this->D_distribution[i].end(), 0);
//    }
//}
//
//double *HeteroRandomConnectivity::GetDistributionJ(long preNeuronId, long postNeuronId) {
//    return &J_distribution[preNeuronId][postNeuronId];
//}
//
//void HeteroRandomConnectivity::SetDistributionJ() {
//    this->J_distribution->resize(this->synapse->GetNoNeuronsPre());
//    for (unsigned long i = 0; i < this->synapse->GetNoNeuronsPre(); ++i) {
//        this->J_distribution[i].resize(this->synapse->GetNoNeuronsPost());
//        std::fill(this->J_distribution[i].begin(), this->J_distribution[i].end(), 1.0);
//    }
//}