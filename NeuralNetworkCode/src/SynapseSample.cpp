#include "SynapseSample.hpp"

SynapseSample::SynapseSample(NeuronPopSample * neurons,std::vector<std::string> *input,GlobalSimInfo *info)
{
    this->generalSynapseSeed = -1;
    this->info    = info;//These are pointers.
    this->neurons = neurons;

    //this->global_D_max = 0;

    int P    = neurons->GetTotalPopulations();

    synapses = new Synapse**[P];
    for(int i = 0; i < P; i++)
        synapses[i] = new Synapse*[P];

    for(int i = 0; i < P; i++){
        for(int j = 0; j < P; j++){
            synapses[i][j] = new CurrentSynapse(neurons->GetPop(i),neurons->GetPop(j),info);
        }
    }

    LoadParameters(input);
}

void SynapseSample::LoadParameters(std::vector<std::string> *input){

    std::string              name,token;
    std::vector<std::string> values;
    int P    = neurons->GetTotalPopulations();

    for(std::vector<std::string>::iterator it = (*input).begin(); it != (*input).end(); ++it) {
        SplitString(&(*it),&name,&values);
        if ((name.find("type") != std::string::npos)&&((name.find("connectivity") == std::string::npos))) {
            SaveSynapseType(name,values.at(0),input);
        }
        if(name.find("generalSynapseSeed") != std::string::npos){
            generalSynapseSeed = std::stoi(values.at(0));
        }
    }

    if(info->globalSeed != -1){
        std::uniform_int_distribution<int> distribution(0,INT32_MAX);
        generalSynapseSeed = distribution(info->globalGenerator);
    }

    if(generalSynapseSeed >= 0){
        //Find P^2 different seed values
        std::mt19937 generator;
        std::uniform_int_distribution<int> distribution(0,INT32_MAX);
        generator = std::mt19937(generalSynapseSeed);

        for(int i = 0; i < P; i++){
            for(int j = 0; j < P; j++){
                synapses[i][j]->SetSeed(&generator);
            }
        }
    }

    //Check seeds (SEPARATE FUNCTION for each pop)
    for(int i = 0; i < P; i++){
        for(int j = 0; j < P; j++){
            for(int k = 0; k < P; k++){
                for(int n = 0; n < P; n++){
                    if((synapses[i][j]->GetSeed() == synapses[k][n]->GetSeed()) && (i>k || j>n)){
                        std::cout << "WARNING! Same Seeds in Connections of Populations " << i << "_" << j << " and " << k << "_" << n << std::endl;
                    }
                }
            }
        }
    }
    /*
    // Get maximal synaptic delay across all synapse populations
    int curr_max = 0;
    for(int i = 0; i < P; i++){
        for(int j = 0; j < P; j++){
            curr_max = synapses[i][j]->GetMaxD();
            if(global_D_max < curr_max){
                global_D_max = curr_max;
            }
        }
    }
    */
}


void SynapseSample::SaveSynapseType(std::string name,std::string type,std::vector<std::string> *input){

    int i,j;
    std::vector<std::string> synapses_strs;
    std::string synNo;
    bool found = false;
    int P      = this->neurons->GetTotalPopulations();

    for(i = 0; i < P; i++){
        for(j = 0; j < P; j++){
            synNo = std::to_string(i) + "_" + std::to_string(j);
            if(name.find(synNo) != std::string::npos && (name.find("connectivity")== std::string::npos)) {
                found = true;
                break;
            }
        }
        if(found)
            break;
    }

    if(!found)
        return;
        //throw("SynapseSample: Synapse id not defined");

    if(type == str_currentSynapse){
        delete synapses[i][j];
        synapses[i][j] = new CurrentSynapse(neurons->GetPop(i),neurons->GetPop(j),info);
    }
    else if(type == str_mongilloSynapse){
        delete synapses[i][j];
        synapses[i][j] = new MongilloSynapse(neurons->GetPop(i),neurons->GetPop(j),info);;
    }
    else if(type == str_mongilloSynapseContinuous){
        delete synapses[i][j];
        synapses[i][j] = new MongilloSynapseContinuous(neurons->GetPop(i),neurons->GetPop(j),info);;
    }
    else if(type == str_prgSynapseContinuous){
        delete synapses[i][j];
        synapses[i][j] = new PRGSynapseContinuous(neurons->GetPop(i),neurons->GetPop(j),info);;
    }
    else if(type == str_exponentialCurrentSynapse) {
        delete synapses[i][j];
        synapses[i][j] = new ExponentialCurrentSynapse(neurons->GetPop(i),neurons->GetPop(j),info);;
    }
	else if (type == str_powerlawsynapse) {
		delete synapses[i][j];
		synapses[i][j] = new PowerLawSynapse(neurons->GetPop(i), neurons->GetPop(j), info);
	}
    else if (type == str_heteroSynapse) {
        delete synapses[i][j];
        synapses[i][j]  = new HeteroCurrentSynapse(neurons->GetPop(i), neurons->GetPop(j), info);
    }

    FilterStringVector(input, "synapses_" + synNo, &synapses_strs);
    synapses[i][j]->LoadParameters(&synapses_strs);
}

void SynapseSample::SaveParameters(std::ofstream * stream){

    std::string syn;
    int P      = this->neurons->GetTotalPopulations();

    *stream <<  "#*************************************************\n";
    *stream <<  "#************** Synaptic Parameters **************\n";
    *stream <<  "#*************************************************\n";
    if(info->globalSeed == -1){
        *stream <<  "synapses_generalSynapseSeed                " << std::to_string(generalSynapseSeed) <<  "\n";
        *stream <<  "#generalSynapseSeed = -1: seeds are defined at individual synapse level.\n";
        *stream <<  "#generalSynapseSeed >= 0: general seed overrides individual seeds.\n";
    }

    for(int i = 0; i < P; i++){
        for(int j = 0; j < P; j++){
            syn = "synapses_" + std::to_string(i) + "_" + std::to_string(j) + "_";
            *stream <<  "#*************************************************\n";
            this->synapses[i][j]->SaveParameters(stream,syn);
        }
    }

}

SynapseSample::~SynapseSample()
{
    int P = this->neurons->GetTotalPopulations();
    for(int i = 0; i < P; i++){
        for(int j = 0; j < P; j++)
        {
            delete synapses[i][j];
        }
        delete [] synapses[i];
    }
    delete [] synapses;
}

void SynapseSample::Test()
{
    int P = this->neurons->GetTotalPopulations();
    for(int source_pop = 0; source_pop < P; source_pop++){
        for(int target_pop = 0; target_pop < P; target_pop++)
            synapses[target_pop][source_pop]->Test();
    }
}

//void SynapseSample::advect(std::vector<std::vector<double>> * synaptic_dV, std::vector<std::vector<std::vector<double>>> * waiting_matrix)

void SynapseSample::advect(std::vector<std::vector<double>> * synaptic_dV){
	int P = this->neurons->GetTotalPopulations();

	for (int source_pop = 0; source_pop < P; source_pop++) {
		for (int target_pop = 0; target_pop < P; target_pop++) {
            if (!synapses[target_pop][source_pop]->IsConnected()){
                continue;
            }
			synapses[target_pop][source_pop]->advect(&synaptic_dV->at(target_pop));
		}

    }
}

void SynapseSample::Reset() {
	int P = this->neurons->GetTotalPopulations();

	for (int source_pop = 0; source_pop < P; source_pop++) {
		for (int target_pop = 0; target_pop < P; target_pop++) {
			synapses[target_pop][source_pop]->ResetWaitingMatrixEntry();
		}

	}
}

void SynapseSample::ConnectNeurons()
{
    int P = this->neurons->GetTotalPopulations();
    for(int sourcePop = 0; sourcePop < P; sourcePop++){
        for(int targetPop = 0; targetPop < P; targetPop++){
            std::cout << "\nConnecting pop. " << sourcePop << " with pop. " << targetPop << ": ";
            this->synapses[targetPop][sourcePop]->ConnectNeurons();
            this->synapses[targetPop][sourcePop]->SetDistributionD();
            this->synapses[targetPop][sourcePop]->SetDistributionJ();
        }
    }

}

void  SynapseSample::WriteConnectivity(std::string filename,int noNeuronsConnectivity){

    std::string filename_ij;
    int P = this->neurons->GetTotalPopulations();
    for(int sourcePop = 0; sourcePop < P; sourcePop++){
        for(int targetPop = 0; targetPop < P; targetPop++){
            filename_ij = filename + "_" + std::to_string(targetPop) + "_" + std::to_string(sourcePop);
            this->synapses[targetPop][sourcePop]->WriteConnectivity(filename_ij,noNeuronsConnectivity);
        }
    }
}

void  SynapseSample::WriteDistributionD(std::string filename,int noNeuronsDelay){

    std::string filename_ij;
    int P = this->neurons->GetTotalPopulations();
    for(int sourcePop = 0; sourcePop < P; sourcePop++){
        for(int targetPop = 0; targetPop < P; targetPop++){
            filename_ij = filename + "_" + std::to_string(targetPop) + "_" + std::to_string(sourcePop);
            this->synapses[targetPop][sourcePop]->WriteDistributionD(filename_ij,noNeuronsDelay);
        }
    }
}

void  SynapseSample::WriteDistributionJ(std::string filename,int noNeuronsJPot){

    std::string filename_ij;
    int P = this->neurons->GetTotalPopulations();
    for(int sourcePop = 0; sourcePop < P; sourcePop++){
        for(int targetPop = 0; targetPop < P; targetPop++){
            filename_ij = filename + "_" + std::to_string(targetPop) + "_" + std::to_string(sourcePop);
            this->synapses[targetPop][sourcePop]->WriteDistributionJ(filename_ij,noNeuronsJPot);
        }
    }
}

int SynapseSample::GetNumberOfDataColumns()
{
    int P            = this->neurons->GetTotalPopulations();
    int data_columns = 0;
    for(int pre = 0; pre < P; pre++){
        for(int post = 0; post <P; post++){
            data_columns += synapses[post][pre]->GetNumberOfDataColumns();
        }
    }
    return data_columns;
}

int SynapseSample::GetNumberOfDataColumns(int post_population, int pre_population)
{
    return synapses[post_population][pre_population]->GetNumberOfDataColumns();
}

std::string SynapseSample::GetDataHeader(int data_column)
{
    std::string data_header;
    int column_counter = 0;
    int              P = this->neurons->GetTotalPopulations();

    for(int post = 0; post < P; post++)
    {
        for(int pre = 0; pre < P ; pre++)
        {
            data_header    += synapses[post][pre]->GetDataHeader(data_column+column_counter);
            column_counter += synapses[post][pre]->GetNumberOfDataColumns();
        }
    }
    return data_header;
}

std::string SynapseSample::GetUnhashedDataHeader()
{
	std::string unhashedheader;
	int              P = this->neurons->GetTotalPopulations();

	for (int post = 0; post < P; post++){
		for (int pre = 0; pre < P; pre++){
			unhashedheader += synapses[post][pre]->GetUnhashedDataHeader();
		}
	}
	return unhashedheader;
}

std::valarray<double> SynapseSample::GetSynapticState(int post_population,int pre_population, int pre_neuron)
{
    return synapses[post_population][pre_population]->GetSynapticState(pre_neuron);
}

double SynapseSample::GetRecurrentInput(int post_population, int pre_population, int post_neuron)
{
	return synapses[post_population][pre_population]->GetrecurrentInput(post_neuron);
}

double SynapseSample::GetCumulatedDV(int post_population, int pre_population)
{
    return synapses[post_population][pre_population]->GetCumulatedDV();
}

int SynapseSample::GetNumberOfPostsynapticTargets(int post_population,int pre_population, int pre_neuron)
{
    return synapses[post_population][pre_population]->GetNumberOfPostsynapticTargets(pre_neuron);
}
