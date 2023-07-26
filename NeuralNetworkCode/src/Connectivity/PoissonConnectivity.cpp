#include "PoissonConnectivity.hpp"
#include "../Synapse/Synapse.hpp"

PoissonConnectivity::PoissonConnectivity(Synapse* synapse, const GlobalSimInfo* infoGlobal):Connectivity(synapse,infoGlobal){
    // std::uniform_int_distribution<int> distribution(0,INT32_MAX);
    // SetSeed(distribution(infoGlobal->globalGenerator));
}

void PoissonConnectivity::SaveParameters(std::ofstream& wParameterStream, std::string idString) const {
    Connectivity::SaveParameters(wParameterStream,idString);
    wParameterStream << idString << "connectivity_connectionLambda\t\t\t" << std::to_string(this->connectionLambda) << "\n";
    // if(infoGlobal->globalSeed == -1){
    //     wParameterStream << idString << "connectivity_seed                  " << std::to_string(this->seed)  << "\n";
    // }
    wParameterStream << "#\t\t"<< IDstringPoissonConnectivity << ": Connectivity between two neurons is determined by a lambda distribution, multiple synapses are enabled.\n";
    // wParameterStream << "#\t\t"<< IDstringBinaryrandomConnectivity << ": Erdoes Renyi network. Each neuronal pair is connected with probability connectionProbability. (as used by [Renart et al. (2010)]). \n";
    //*stream << "# "<< stringBinaryrandomConnectivity << "as used by [Renart et al. (2010)] \n";
}

void PoissonConnectivity::LoadParameters(const std::vector<FileEntry>& connectivityParameters){
    Connectivity::LoadParameters(connectivityParameters);

    for(auto& [parameterName, parameterValues] : connectivityParameters) {
        // if(parameterName.find("seed") != std::string::npos){
        //     SetSeed(std::stoi(parameterValues.at(0)));
        // }
        if (parameterName.find("ConnectionProba") != std::string::npos || parameterName.find("connectionLambda") != std::string::npos) {
            this->connectionLambda = std::stod(parameterValues.at(0));
        }
    }
}

void PoissonConnectivity::ConnectNeurons()
{
    std::poisson_distribution<> poissonDistribution(connectionLambda);
    signed long  noConnectionsFormed;
    NeuronInt    noTargetNeurons = synapse->GetNoTargetNeurons();
	NeuronInt    noSourceNeurons = synapse->GetNoSourceNeurons();

    NeuronInt    outputInterval = noTargetNeurons/10;
    if(outputInterval == 0){
        outputInterval = 1;
    }

    //std::uniform_real_distribution<double> distribution(0,1);


    //Iterate through all target neurons
    for(NeuronInt targetNeuron : std::ranges::views::iota(0,noTargetNeurons)){

        //Iterate through all source neurons
        for(NeuronInt sourceNeuron : std::ranges::views::iota(0,noSourceNeurons)){

            //Connect with given probability
            noConnectionsFormed = poissonDistribution(generator);
            totalConnections += noConnectionsFormed;
            if (noConnectionsFormed==0){
                continue;
            }
            if((synapse->IsRecurrent()) && (sourceNeuron == targetNeuron)){
                continue;
            }
            //Check if target was assigned to the same sourceNeuron already
            // if(!synapse->IsSourceVectorEmpty(sourceNeuron) && synapse->IsTargetLastInVector(targetNeuron,sourceNeuron)){
            //     continue;
            // }

            for (int connection : std::ranges::views::iota(0, noConnectionsFormed)){
                (void)connection;
                synapse->AllocateSynapse(targetNeuron, sourceNeuron);
            }
        }

        // if(targetNeuron%outputInterval == 0)
        //     std::cout << 100*(targetNeuron)/noTargetNeurons << "%\n";
    }
        // std::cout << "100%\n";
}


double PoissonConnectivity::GetExpectedConnections() const{
    // return connectionLambda*this->synapse->GetNoSourceNeurons();
    return totalConnections/synapse->GetNoTargetNeurons();
}
