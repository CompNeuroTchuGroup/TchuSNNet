#include "NeuronPop.hpp"

NeuronPop::NeuronPop(GlobalSimInfo * info,int id){
    refractorySteps        = 0;
    this->info             = info;
    seed_InitialPotentials = 0;
    seed_InitialPrevSpike  = 0;
    noNeurons              = 0;
    this->identifier       = id;

    tau_m                  = 0;
    v_reset                = 0;
    v_thresh               = 0;
}


void NeuronPop::SetSeeds(int seed1,int seed2) {
    if(info->globalSeed != -1){
        std::uniform_int_distribution<int> distribution(0,INT32_MAX);
        seed_InitialPrevSpike = distribution(info->globalGenerator);
        seed_InitialPotentials = distribution(info->globalGenerator);
    } else {
        seed_InitialPrevSpike  = seed1;
        seed_InitialPotentials = seed2;
    }
}


void NeuronPop::SetNeurons(unsigned long noNeur) {
    noNeurons = noNeur;

    srand(seed_InitialPrevSpike);

    previous_spike_step.resize(noNeurons);
    potential.resize(noNeurons);

    //std::cout << "v_thresh = " << std::to_string(this->v_thresh) << "\n";
    //std::cout << "v_reset = " << std::to_string(this->v_reset) << "\n";
    
    std::mt19937 generator(seed_InitialPotentials);
    //std::uniform_real_distribution<double> uniformDistribution (0.0,1.0); //v_thresh instead of 1.0
    std::uniform_real_distribution<double> uniformDistribution (v_reset,v_thresh);

    for(unsigned long i = 0; i < noNeurons; i++){
        previous_spike_step[i] = - rand() % static_cast<int>(3.0/info->dt);
        potential[i]           = uniformDistribution(generator);
    }
}

void NeuronPop::SetPosition(long noNeur)
{
	if (info->Dimensions == 0)
		return;
	noNeurons = noNeur;
	std::mt19937 generator(seed_InitialPotentials);
	std::uniform_real_distribution<double> uniformDistribution(0.0, 1.0);
	x_pos.resize(noNeurons);
	y_pos.resize(noNeurons);
	double Ly = info->Ly;
	double Lx = info->Lx;
	int rows;
	int columns;
	int mod;

	//for (int i = 0; i < noNeurons; i++) {
	//	x_pos[i] = uniformDistribution(generator) * info->Lx;
	//	y_pos[i] = uniformDistribution(generator) * info->Ly;
	//}

	if (info->Dimensions == 2) {
        rows = static_cast<int>(ceil(sqrt(noNeurons)));
		columns = noNeurons / rows;
		mod = noNeurons -rows*columns;//the difference is added by including one extra column in each of the first mod rows

		for (int i = 0; i < mod*(columns+1); i++) {
			x_pos[i] = (Lx / (columns + 1))*(i % (columns + 1));
			y_pos[i] = (Ly / rows) *floor(i / (columns + 1));
		}
		for (unsigned long i = (mod+mod*columns); i < noNeurons; i++) {
			x_pos[i] = (Lx / (columns))*((i - mod) % columns);
			y_pos[i] = (Ly / rows) *floor((i - mod) / (columns));
		}
	}
	else if (info->Dimensions == 1) {
		for (unsigned long i = 0;i < noNeurons; i++) {
			x_pos[i] = i * Lx / noNeurons;
		}
	}
}

//void NeuronPop::SetGlobalNeuronId(long global_neuron_id){
//    global_id.resize(noNeurons);
//    for(long i = 0; i < noNeurons; i++){
//        global_id[i] = global_neuron_id + i;
//    }
//}

void NeuronPop::ClearSpiker(){
    for(const auto &previous_spiker : spiker)
        previous_spike_step[previous_spiker] = info->time_step-1;

    spiker.clear();
}

void NeuronPop::LoadParameters(std::vector<std::string> *input,unsigned long noNeurons){
    SetNeurons(noNeurons);
    LoadParameters(input);
}

void NeuronPop::LoadParameters(std::vector<std::string> *input){

    std::string              name,token;
    std::vector<std::string> values;

    for(std::vector<std::string>::iterator it = (*input).begin(); it != (*input).end(); ++it) {

        SplitString(&(*it),&name,&values);

        /*if(name.find("Ni") != std::string::npos){
            totalPopulations = (int)values.size();
            //spiker              = new std::vector<int>[totalPopulations];
            //neuronsInPopulation = new int[totalPopulations];
            //for(int i = 0;i<totalPopulations;i++)
            //    neuronsInPopulation[i] = std::stoi(values.at(i));
        }*/
        if(name.find("tauM") != std::string::npos){
            this->tau_m = (std::stod(values.at(0)));
        }
        else if(name.find("vReset") != std::string::npos){
            this->v_reset = (std::stod(values.at(0)));
        }
        else if(name.find("vThresh") != std::string::npos){
            this->v_thresh = (std::stod(values.at(0)));
        }
        else if(name.find("refractoryTime") != std::string::npos){
            this->refractorySteps = std::round((std::stod(values.at(0)))/info->dt);
        }
        //else if(name.find("r_target") != std::string::npos){
         //   this->r_target = (std::stoi(values.at(0)));
        //}
        else if(name.find("seedInitialPotentials") != std::string::npos){
            this->seed_InitialPotentials = (std::stoi(values.at(0)));
        }
        else if(name.find("seedInitialPrevSpike") != std::string::npos){
            this->seed_InitialPrevSpike = (std::stoi(values.at(0)));
        }
        else if(name.find("noNeurons") != std::string::npos){
                    SetNeurons(static_cast<long>(std::stod(values.at(0))));
        }
        // else if(name.find("streamOutput") != std::string::npos){
        //             if (values.at(0).find("true") != std::string::npos){
        //                 streamingNOutputBool=true;
        //             }
        // }
    }

    // for(std::vector<std::string>::iterator it = (*input).begin(); it != (*input).end(); ++it) {
    //
    //     SplitString(&(*it),&name,&values);
    //
    //     if(name.find("noNeurons") != std::string::npos){
    //         SetNeurons((long)(std::stod(values.at(0))));
    //         //spiker              = new std::vector<int>[totalPopulations];
    //         //neuronsInPopulation = new int[totalPopulations];
    //         //for(int i = 0;i<totalPopulations;i++)
    //         //    neuronsInPopulation[i] = std::stoi(values.at(i));
    //     }
    // }
}

void NeuronPop::SaveParameters(std::ofstream * stream){

    std::string id = "neurons_" + std::to_string(GetId());

    *stream <<  "#***********************************************\n";
    //*stream <<  id + "_streamOutput                " << std::boolalpha<< streamingNOutputBool << std::noboolalpha << "\n";
    *stream <<  id + "_noNeurons                   " << noNeurons << "\n";
    *stream <<  id + "_type                        " << GetType() << "\n";
    *stream <<  id + "_tauM                        " << std::to_string(this->tau_m)  << " #seconds\n";
    *stream <<  id + "_vReset                      " << std::to_string(this->v_reset)  << " #mV \n";
    *stream <<  id + "_vThresh                     " << std::to_string(this->v_thresh)  << " #mV\n";
    *stream <<  id + "_refractoryTime              " << std::to_string(this->refractorySteps*info->dt)  << " #seconds\n";
    //*stream <<  id + "_r_target                    " << std::to_string(this->r_target)  << " Hz\n";
    if(info->globalSeed == -1){
        *stream <<  id + "_seedInitialPotentials   " << this->seed_InitialPotentials << "\n";
        *stream <<  id + "_seedInitialPrevSpike    " << this->seed_InitialPrevSpike << "\n";
    }
    *stream <<  "#\t\tNote: Resting potential is 0 by definition.\n";


}


//From here on 

BaseSpinePtr NeuronPop::AllocateNewSynapse(unsigned long neuronId, HeteroCurrentSynapse&syn) {
        assertm(false, "Non-hetero NeuronPop called AllocateNewSynapse");
        throw; 
        BaseSpinePtr empty{};
        return empty;
}

void NeuronPop::RecordExcitatorySynapticSpike(unsigned long neuronId, unsigned long synapseId){
    throw;
}

std::valarray<double> NeuronPop::GetOverallSynapticProfile(unsigned long neuronId){
    throw;
    std::valarray<double> empty{};
    return empty;
}

std::valarray<double> NeuronPop::GetIndividualSynapticProfile(unsigned long neuronId, unsigned long synapseId){
    throw;
    std::valarray<double> empty{};
    return empty;
}