
#include "Synapse.hpp"

Synapse::Synapse(NeuronPop * postNeurons, NeuronPop * preNeurons,GlobalSimInfo * info){
    this->info         = info;
    this->neuronsPre   = preNeurons;
    this->neuronsPost  = postNeurons;
    //neurons            = ptr;
    //pre_population_id  = preId;
    //post_population_id = postId;
    geometry           = new RandomConnectivity(this,info);
    cumulatedDV        = 0.0;
    //********************************
    //***** Set Default values *******
    //********************************
    J                  = 0;
    SigmaJ             = 0;
    J_pot              = 0;
    P_pot              = 0;

    D_min  = 0;
    D_max = 0;
    //********************************
    //********************************
}

const std::string Synapse::GetIdStr(){
    int preId  = this->neuronsPre->GetId();
    int postId = this->neuronsPost->GetId();
    std::string idStr = std::to_string(postId) + std::to_string(preId);
    return idStr;
}

std::string Synapse::GetIdStrWithULine() {
    int preId = this->neuronsPre->GetId();
    int postId = this->neuronsPost->GetId();
    std::string idStr = std::to_string(postId) + "_" + std::to_string(preId);
    return idStr;
}

void Synapse::LoadParameters(std::vector<std::string> *input){

    std::vector<std::string> connectivity_strs;
    std::string              name;
    std::vector<std::string> values;
    double                   dt = info->dt;
    //std::cout << "loading parameters for synapse " << std::to_string(this->neuronsPre->GetId()) <<"->" <<std::to_string(this->neuronsPost->GetId()) << "\n";

    // These boolean variables keeps track of whether J_pot and P_pot are specified in the parameter file. If not, they are set to J and 0 respectively after the loop.
    bool found_J_pot = false;
    bool found_P_pot = false;

    for(std::vector<std::string>::iterator it = (*input).begin(); it != (*input).end(); ++it) {
        SplitString(&(*it),&name,&values);
        if(name.find("J_pot") != std::string::npos){
            this->J_pot = std::stod(values.at(0));
            found_J_pot = true;
            //std::cout << "assigning J_pot = " << std::to_string(this->J_pot) << "\n";
        }
        else if(name.find('J') != std::string::npos){
            this->J = std::stod(values.at(0));
            //std::cout << "assigning J = " << std::to_string(this->J) << "\n";
        }
        else if (name.find("Sigma_j") != std::string::npos) {
            this->SigmaJ = std::stod(values.at(0));
            //std::cout << "assigning SigmaJ = " << std::to_string(this->SigmaJ) << "\n";
        }
        else if(name.find("P_pot") != std::string::npos){
            this->P_pot = std::stod(values.at(0));
            found_P_pot = true;
            //std::cout << "assigning P_pot = " << std::to_string(this->P_pot) << "\n";
        }
        else if(name.find("D_min") != std::string::npos){
            this->D_min = static_cast<int>(std::round(std::stod(values.at(0)) / dt));
        }
        else if(name.find("D_max") != std::string::npos){
            this->D_max = static_cast<int>(std::round(std::stod(values.at(0)) / dt));
        }
        else if(name.find("connectivity_type") != std::string::npos){
            if(values.at(0) == str_randomConnectivity){
                if(geometry != nullptr) //Doing NOTHING, delete already checks
                    delete geometry;
                geometry           = new RandomConnectivity(this,info);
            }
            else if(values.at(0) == str_binaryrandomConnectivity){
                if(geometry != nullptr)
                    delete geometry;
                geometry       = new BinaryRandomConnectivity(this,info);
            }
			else if (values.at(0) == str_distanceConnectivity) {
				if (geometry != nullptr)
					delete geometry;
				geometry = new DistanceConnectivity(this, info);
			} else if(values.at(0) == str_heteroRandomConnectivity) {
			    delete geometry;
                this->geometry = new HeteroRandomConnectivity(this, info);
            } else if (values.at(0) == str_adjacencyMatrixConnectivity) {
                delete geometry;
                this->geometry = new AdjacencyMatrixConnectivity(this, info);
            }
        } else if(name.find("connected") != std::string::npos){
            if (values.at(0).find("false")!=std::string::npos){
                this->SetConnectedFalse();
            } else if (values.at(0).find("true")!=std::string::npos){
                this->isConnectedBool=true; //just to make sure
            } else {
                throw; //if connected is there, but the input is not correct, exception
            }
        }
    }

    // Set default values of the optional parameters J_pot and P_pot if necessary

    if (!found_P_pot){
        this->P_pot = 0;
    }
	else if (!found_J_pot) {
		this->J_pot = this->J;
	}

	/*
	std::vector<std::vector<std::vector<double>>> waiting_matrix (global_D_max+1, std::vector<std::vector<double>>(P, std::vector<double>(0)));
	for(int d_index = 0; d_index < waiting_matrix.size(); d_index++){
	for(int p_index = 0; p_index < waiting_matrix[0].size(); p_index++){
	waiting_matrix[d_index][p_index].resize(neurons->GetNeuronsPop(p_index),0);
	}
	}*/

	waiting_matrix.resize(GetNoNeuronsPost());
	for (unsigned long i = 0;i < GetNoNeuronsPost();i++)
		waiting_matrix[i].resize(static_cast<size_t>(D_max) + 1);

    FilterStringVector(input,"connectivity", &connectivity_strs); //to use once all input Parameter files are updated
    geometry->LoadParameters(&connectivity_strs);
}

void Synapse::SaveParameters(std::ofstream * stream,std::string id_str){
    *stream << id_str << "type\t\t\t\t\t\t\t" << GetTypeStr() << "\n";
    *stream << id_str << "connected\t\t\t\t\t\t" << std::boolalpha << this->isConnectedBool <<std::noboolalpha<< "\n";
    *stream << id_str << "D_min\t\t\t\t\t\t" << std::to_string(this->D_min*info->dt) << " seconds\n";
    *stream << id_str << "D_max\t\t\t\t\t\t" << std::to_string(this->D_max*info->dt) << " seconds\n";
    *stream << id_str << "J\t\t\t\t\t\t\t" << std::to_string(this->J) << " dmV/Spike\n";
    *stream << id_str << "Sigma_j\t\t\t\t\t\t" << std::to_string(this->SigmaJ) << " dmV/Spike\n";
    *stream << id_str << "J_pot\t\t\t\t\t\t" << std::to_string(this->J_pot) << " dmV/Spike\n";
    *stream << id_str << "P_pot\t\t\t\t\t\t" << std::to_string(this->P_pot) << "\n";

    geometry->SaveParameters(stream,id_str);
}


unsigned long Synapse::GetNumberOfPostsynapticTargets(int pre_neuron) {
    return static_cast<unsigned long>(geometry->GetTargetList(pre_neuron)->size());
}

double Synapse::GetrecurrentInput(int post_neuron){
	int index = (this->info->time_step) % (D_max + 1);
	return waiting_matrix[post_neuron][index];
}

void Synapse::advect(std::vector<double> * synaptic_dV)
{
    resetcumulatedDV();//Does this do something ?
    std::vector<long>  *    spikersUsed;
    std::vector<double>     currents;

    //Get list of spikers
    spikersUsed = neuronsPre->GetSpikers();

    //Go through all the spikers and add current arising from spikers to waiting_matrix
    for(auto const& spiker: *spikersUsed){
        std::vector<unsigned long> tL = *geometry->GetTargetList(spiker);
        currents.resize(tL.size(),0);
        std::fill(currents.begin(), currents.end(), 0);//OPTIMIZATION: Change this to the HCS function (currents)

        advectSpikers(currents, spiker);//Here target list is used again OPTIMIZATION. Both these optimizations can be done, but it requires a revision of all Synapse classes.


        FillWaitingMatrix(spiker, currents);//Here target list is used again OPTIMIZATION
    }

	ReadWaitingMatrixEntry(*synaptic_dV);
	//advect_finalize(&waiting_matrix);// Is it OK for Conductance, Exponential and Graupner Synapse ? Fill Waiting Matrix has already been called
}

void Synapse::FillWaitingMatrix(long spiker, std::vector<double>& currents){
    std::vector<unsigned long> tL {*geometry->GetTargetList(spiker)};
    //std::vector<int> *tD = geometry->GetDistributionD(spiker);

    for(unsigned int target{0}; target < tL.size(); ++target){
        int delay { *geometry->GetDistributionD(spiker,target)};                              // get individual synaptic delay between spiker and target
        int matrix_index { (this->info->time_step + delay)%(D_max+1)};
        waiting_matrix[tL[target]][matrix_index] += currents[target]; // add spiker to waiting matrix
    }
}

void Synapse::resetcumulatedDV(){
	cumulatedDV = 0.0;
}

void Synapse::ResetWaitingMatrixEntry() {
	int index= (this->info->time_step) % (D_max + 1);
	for (unsigned long i = 0;i < GetNoNeuronsPost();i++)
		(waiting_matrix)[i][index]=0;
}

void Synapse::ReadWaitingMatrixEntry(std::vector<double>& synaptic_dV) {
	int index = (this->info->time_step) % (D_max + 1);
	for (unsigned long i = 0;i < GetNoNeuronsPost();i++)
		synaptic_dV[i] = synaptic_dV[i]+ waiting_matrix[i][index];
}

void Synapse::SetSeed(std::mt19937 *generator){

    geometry->SetSeed(generator);
    //std::uniform_int_distribution<int> distribution(0,INT32_MAX);
    //std::cout << distribution(*generator) << "\n";
}

void Synapse::ConnectNeurons(){
    if (IsConnected()){
        geometry->ConnectNeurons();
    } else {
        std::cout<< "Connection skipped." <<"\n";
    }
    localScalingFactor=pow(geometry->GetNumberAverageSourceNeurons(),info->networkScaling_synStrength);
    globalScalingFactor=pow(info->N,info->networkScaling_synStrength);
}

/*double Synapse::GetCouplingStrength(){
    if(info->networkScaling_extern == 0){
        return (J*pow(geometry->GetNumberAverageSourceNeurons(),info->networkScaling_synStrength));}
    else if(info->networkScaling_extern == 1)
        return (J*pow(info->N,info->networkScaling_synStrength));
    else
        throw "error in Synapse::GetCouplingStrength";
}*/

double Synapse::GetCouplingStrength(unsigned long preNeuronId, unsigned long postNeuronId){ // For distribution of Js
    J = *geometry->GetDistributionJ(preNeuronId,postNeuronId);
    if(info->networkScaling_extern == 0){
        return (J*localScalingFactor);}
    else if(info->networkScaling_extern == 1)
        {return (J*globalScalingFactor);}
    else
        {throw "error in Synapse::GetCouplingStrength";}
}

void Synapse::SetDistributionD(){
    geometry->SetDistributionD();
}

void Synapse::SetDistributionJ(){
    geometry->SetDistributionJ();
}

NeuronPop* Synapse::GetNeuronsPre() {
    return neuronsPre;
}
// std::vector<std::vector<double>> getWaitingMatrix(const Synapse& synapse) {
//     std::vector<std::vector<double>> copy(synapse.waiting_matrix.size());
//     for (int i = 0; i < synapse.waiting_matrix.size(); ++i) {
//         for (auto& item: synapse.waiting_matrix[i]) {
//             copy[i].push_back(item);
//         }
//     }
//     return copy;
// }

NeuronPop* Synapse::GetNeuronsPost() {
    return neuronsPost;
}