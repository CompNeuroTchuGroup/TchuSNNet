#include "RandomConnectivity.hpp"
#include "../Synapse/Synapse.hpp"

RandomConnectivity::RandomConnectivity(Synapse * syn,GlobalSimInfo  * info):Connectivity(syn,info){
    noSourceNeurons    = 0;
    SetSeed(0);
}

void RandomConnectivity::SaveParameters(std::ofstream * stream,std::string id_str){
    Connectivity::SaveParameters(stream,id_str);
    //*stream << id_str << "connectivity_noSourceNeurons " << std::to_string(this->noSourceNeurons) << "\n";
    *stream << id_str << "connectivity_ConnectionProba\t\t\t" << std::to_string(this->GetConnectionProbability()) << "\n";
    //*stream << "#" << id_str << "connectivity_noSourceNeurons " << std::to_string(this->noSourceNeurons) << "\n";
    *stream << "#\t\t"<< str_randomConnectivity << ": Each neuron receives C = connectionProbability*N_p randomly chosen connections from the presynaptic population p (as used by [Brunel (2000)]).\n";
}

void RandomConnectivity::LoadParameters(std::vector<std::string> *input){
    Connectivity::LoadParameters(input);

    std::string              name;
    std::vector<std::string> values;

    for(std::vector<std::string>::iterator it = (*input).begin(); it != (*input).end(); ++it) {
        SplitString(&(*it),&name,&values);
        if(name.find("seed") != std::string::npos){
            SetSeed(std::stoi(values.at(0)));
        }
		else if (name.find("ConnectionProba") != std::string::npos || name.find("connectionProbability") != std::string::npos) {
            //std::cout << "number of neurons pre: " << std::to_string(synapse->GetNoNeuronsPre()) << ", ";
            //std::cout << "connection probability value: " << std::to_string(std::stod(values.at(0))) << ", ";
            //std::cout << "their product: " << std::to_string(synapse->GetNoNeuronsPre()*std::stod(values.at(0))) << "\n";
            SetNoSourceNeurons(static_cast<unsigned long>(std::lround(synapse->GetNoNeuronsPre() * std::stod(values.at(0)))));
        }
        else if(name.find("noSourceNeurons") != std::string::npos){
            SetNoSourceNeurons(std::stoi(values.at(0)));
        }

    }
}

void RandomConnectivity::SetNoSourceNeurons(unsigned long noSN){
    if(noSN < 0)
        noSourceNeurons = 0;
    else if(noSN > synapse->GetNoNeuronsPre())
        noSourceNeurons = synapse->GetNoNeuronsPre();
    else if( (noSN == synapse->GetNoNeuronsPre()) && (synapse->IsRecurrent()))
        //noSourceNeurons = 0;
        noSourceNeurons = synapse->GetNoNeuronsPre()-1;
    else
        noSourceNeurons = noSN;

}

double RandomConnectivity::GetConnectionProbability(){
    if(synapse ->GetNoNeuronsPre() == 0)
        return 0;
    if((noSourceNeurons == synapse->GetNoNeuronsPre()) && synapse->IsRecurrent()){
        return (static_cast<double>(noSourceNeurons+1))/(static_cast<double>(synapse ->GetNoNeuronsPre()));
    }
    else
        return (static_cast<double>(noSourceNeurons))/(static_cast<double>(synapse ->GetNoNeuronsPre()));
}


void RandomConnectivity::ConnectNeurons()
{

    unsigned long    source,countedSourceNeurons;
    unsigned long    nPost = synapse->GetNoNeuronsPost();

    long    output_Interval = nPost/10;
    if(output_Interval == 0)
        output_Interval = 1;

    std::uniform_int_distribution<int> distribution(0,(int)synapse->GetNoNeuronsPre()-1);

    //every target neuron has a fixed number of source neurons

    //Iterate through all target neurons
    for(unsigned long target = 0; target < nPost; target++){

        countedSourceNeurons = 0;

        //assign for each target neuron 'noSourceNeurons' source neurons
        while(countedSourceNeurons < noSourceNeurons)
        {
            source = distribution(generator);

            //Check if source and target neurons are equal
            if((synapse->IsRecurrent()) && (source == target))
                continue;

            //Check if target was assigned to the same source already
            if(!target_id[source].empty() && target_id[source].back() == target)
                continue;

            target_id[source].push_back(target);
            countedSourceNeurons++;
        }
        if((target)%output_Interval == 0)
            std::cout << 100*(target)/nPost << "%-";
    }
    std::cout << "100%\n";
}
