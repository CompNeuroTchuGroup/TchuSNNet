
#include "Connectivity.hpp"
#include "../Synapse/Synapse.hpp"

Connectivity::Connectivity(Synapse * syn,GlobalSimInfo  * info){
    this->info = info;
    synapse   = syn;

    // SUGGESTION: Use vector of vectors instead of array of vectors
    // https://github.com/saiftyfirst/BP_Demos/blob/master/C%2B%2B/arrayOfVectors_vs_vectorOfVectors.cpp
    target_id = new std::vector<unsigned long> [synapse->GetNoNeuronsPre()];
    D_distribution = new std::vector<int> [synapse->GetNoNeuronsPre()];
    J_distribution = new std::vector<double> [synapse->GetNoNeuronsPre()];

//	HasPot = true;
    HasJDistribution = true; //Either (Jpot and Ppot) or Sigma If not HasJDistribution, then all synapses have the same j
	HasDdistribution = true;//By default, set as true

    fixSeed = false; //by default, seed can be changed
}

Connectivity::~Connectivity(){
    delete [] target_id;
    delete [] D_distribution;
    delete [] J_distribution;
}

void Connectivity::SaveParameters(std::ofstream* stream, std::string id_str) {
    *stream << id_str << "connectivity_type\t\t\t\t" << GetTypeStr() << "\n";
    //if(info->globalSeed == -1){
    *stream << id_str << "connectivity_seed\t\t\t\t" << std::to_string(this->seed)  << "\n";
    //}
}

void Connectivity::Test(){

    for(unsigned long i = 0;i<synapse->GetNoNeuronsPre();i++){
        for(unsigned long j : target_id[i]){
            //std::cout << synapse->GetNoNeuronsPost() << " -- ";
            if(j >= synapse->GetNoNeuronsPost()){
                std::cout << synapse->GetIdStr() << " - ";
                std::cout << std::to_string(j) << "\n";
                return;
            }
        }
    }
    std::cout << synapse->GetIdStr() << " - Test OK \n";
}

void  Connectivity::SetSeed(std::mt19937 *seedGenerator)  {
    if(!fixSeed){
        std::uniform_int_distribution<int> distribution(0,INT32_MAX);
        SetSeed(distribution(*seedGenerator));
    }
}

void  Connectivity::SetSeed(int s){
    if(!fixSeed){
        seed       = s;
        generator  = std::mt19937(seed);
    }
}

void Connectivity::LoadParameters(std::vector<std::string> *input){
    std::string              name;
    std::vector<std::string> values;

    for(auto & it : *input) {
        SplitString(&it,&name,&values);
        if(name.find("seed") != std::string::npos){
            SetSeed(std::stoi(values.at(0)));
            fixSeed = true;
        }
    }

}

void Connectivity::WriteConnectivity(const std::string& filename,unsigned long noNeuronsConnectivity){

    int count;
    unsigned long i;
    unsigned long noTargetNeurons = synapse->GetNoNeuronsPost();
    unsigned long noSourceNeurons = synapse->GetNoNeuronsPre();

    if(noNeuronsConnectivity <= 0)
        return;

    std::ofstream stream(filename + ".txt"); //ios::binary

    if(noTargetNeurons > noNeuronsConnectivity)
        noTargetNeurons = noNeuronsConnectivity;

    if(noSourceNeurons > noNeuronsConnectivity)
        noSourceNeurons = noNeuronsConnectivity;

    //iterate through source & target neurons & write connectivity
    for(unsigned long source = 0;source < noSourceNeurons;source++){

        i = 0;
        for(unsigned long target = 0;target<noTargetNeurons;target++){
			if (target_id[source].size() <= i || target_id[source][i] != target)
				stream << "0 ";
			else {
				count = 1;
				i++;
				while (target_id[source].size() > i && target_id[source][i] == target) { //in case the same pair can be connected multiple times
					i++;
					count++;
				}
				stream << std::to_string(count)<<" ";
			}
        }
        stream << "\n";

        /*stream << "Targets for neuron " << source << ":" << target_id[source].size() << " -- ";
        for(auto const& target: (target_id[source])){
            stream << std::to_string(target) << " ";
        }
        stream << "\n";*/
    }
    stream.close();

}

int * Connectivity::GetDistributionD(long preNeuronId, long postNeuronId){
	if (HasDdistribution)
		return &D_distribution[preNeuronId][postNeuronId];
	else {
		return (synapse->GetDpointer());
	}

}

void Connectivity::SetDistributionD(){
    std::uniform_int_distribution<int> uniformDistribution(this->synapse->GetMinD(),this->synapse->GetMaxD());
	unsigned long noTarget;

	if (synapse->GetMaxD()== synapse->GetMinD()) {
		HasDdistribution = false;
		//return;
	}
    //std::cout << "printing D_distribution: \n";
    for(unsigned long source_neuron = 0; source_neuron < this->synapse->GetNoNeuronsPre(); source_neuron++){
        if(target_id[source_neuron].size()!=0){
            noTarget = static_cast<unsigned long>((target_id[source_neuron]).size());
			D_distribution[source_neuron].resize(noTarget);
			for (unsigned long target_neuron = 0; target_neuron < noTarget; target_neuron++) {
                int d = uniformDistribution(generator);
                D_distribution[source_neuron][target_neuron]=d;
                //std::cout << std::to_string(d) << " ";
            }
            //std::cout << "\n";
        }
    }
}

double * Connectivity::GetDistributionJ(long preNeuronId, long postNeuronId)
{
	if (HasJDistribution)
		return &J_distribution[preNeuronId][postNeuronId];
	else {
		return (synapse->GetJpointer());
	}

}

void Connectivity::SetDistributionJ(){
    std::uniform_real_distribution<double> real_distribution(0,1);
    std::normal_distribution<double> Gaussian(0, 1);
	unsigned long noTarget;
    double rdj;//random number to determine if the synapse is potentiated
    double rds;//random number to determine the deviation to the mean j
    double Threshold;
    double j;


	if (synapse->GetSigmaJ() == 0 && (synapse->GetPPot() == 0 || this->synapse->GetJPot() == this->synapse->GetJBase())) {
        HasJDistribution = false;
		return;
	}
    //std::cout << "printing J_distribution: \n";
    Threshold = this->synapse->GetPPot();
    for(unsigned long source_neuron = 0; source_neuron < this->synapse->GetNoNeuronsPre(); source_neuron++){
        noTarget = static_cast<unsigned long>((target_id[source_neuron]).size());
		J_distribution[source_neuron].resize(noTarget);
		//std::cout <<  "\n";
        for(unsigned long target_neuron = 0; target_neuron < noTarget; target_neuron++){
			rdj = real_distribution(generator);
            rds = Gaussian(generator);
            if(rdj < Threshold){
				j = (this->synapse->GetJPot());
            }
            else{
				j = (this->synapse->GetJBase());
            }
            j = j + synapse->GetSigmaJ() * rds;
            //if ((synapse->GetJPot() * synapse->GetJBase() >= 0) && (j * synapse->GetJBase() < 0)) {
            //    j = 0;//Dale's law will be strictly enforced unless the provided J and Jpot already blatantly violate it
            //}
            J_distribution[source_neuron][target_neuron] = j;
		}
		//std::cout << "\n";
    }
}


void Connectivity::WriteDistributionD(const std::string& filename, unsigned long noNeuronsDelay){

    int count;
    unsigned long i;
    double delay;
    unsigned long noTargetNeurons = synapse->GetNoNeuronsPost();
    unsigned long noSourceNeurons = synapse->GetNoNeuronsPre();

    if(noNeuronsDelay <= 0)
        return;

    std::ofstream stream(filename + ".txt"); //ios::binary

    if(noTargetNeurons > noNeuronsDelay)
        noTargetNeurons = noNeuronsDelay;

    if(noSourceNeurons > noNeuronsDelay)
        noSourceNeurons = noNeuronsDelay;

    //iterate through source & target neurons & write connectivity
    for(unsigned long source = 0;source < noSourceNeurons;source++){

        i = 0; // counts the number of connections that have been written for this source neuron
        for(unsigned long target = 0;target<noTargetNeurons;target++){
			if (target_id[source].size() <= i || target_id[source][i] != target) // if all connections have been written OR target id does not match
				stream << "nan ";
			else {
				count = 1;
				if (HasDdistribution) {
					delay = D_distribution[source][i];
				}
				else {
					delay = synapse->GetMaxD();
				}
				i++;
				while (target_id[source].size() > i && target_id[source][i] == target) { // while connections are still to be written and target id matches
					i++;
					count++;
					if (HasDdistribution) {
						delay += D_distribution[source][i];
					}
					else {
						delay += synapse->GetMaxD();
					}
				}
				stream << std::to_string(delay/count)<<" ";
			}
        }
        stream << "\n";

        /*stream << "Targets for neuron " << source << ":" << target_id[source].size() << " -- ";
        for(auto const& target: (target_id[source])){
            stream << std::to_string(target) << " ";
        }
        stream << "\n";*/
    }
    stream.close();

}


void Connectivity::WriteDistributionJ(const std::string& filename,unsigned long noNeuronsJPot){

    int count;
    unsigned long i;
    double J;
    unsigned long noTargetNeurons = synapse->GetNoNeuronsPost();
    unsigned long noSourceNeurons = synapse->GetNoNeuronsPre();

    if(noNeuronsJPot <= 0)
        return;

    std::ofstream stream(filename + ".txt"); //ios::binary

    if(noTargetNeurons > noNeuronsJPot)
        noTargetNeurons = noNeuronsJPot;

    if(noSourceNeurons > noNeuronsJPot)
        noSourceNeurons = noNeuronsJPot;

    //iterate through source & target neurons & write connectivity
    for(unsigned long source = 0;source < noSourceNeurons;source++){
        i = 0; // counts the number of connections that have been written for this source neuron
        for(unsigned long target = 0;target<noTargetNeurons;target++){
			if (target_id[source].size() <= i || target_id[source][i] != target) // if all connections have been written OR target id does not match
				stream << "nan ";
			else {
				count = 1;
				if (HasJDistribution) {
					J = J_distribution[source][i];
				}
				else {
					J = synapse->GetJBase();
				}
				i++;
				while (target_id[source].size() > i && target_id[source][i] == target) { // while connections are still to be written and target id matches
					i++;
					count++;
					if (HasJDistribution) {
						J += J_distribution[source][i];
					}
					else {
						J += synapse->GetJBase();
					}
				}
				stream << std::to_string(J/count)<<" ";
			}
        }
        stream << "\n";

        /*stream << "Targets for neuron " << source << ":" << target_id[source].size() << " -- ";
        for(auto const& target: (target_id[source])){
            stream << std::to_string(target) << " ";
        }
        stream << "\n";*/
    }
    stream.close();

}

//After this, all functions are virtual, not intended for running
std::vector<std::pair<unsigned long, unsigned long>>& Connectivity::getSynapticTargets(const unsigned long&){
    throw;
    std::vector<std::pair<unsigned long, unsigned long>> empty{};
    return empty;
    }