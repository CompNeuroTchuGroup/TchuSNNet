#include "NeuralNetwork.hpp"

NeuralNetwork::NeuralNetwork(std::string baseDir,std::vector<FileEntry> *parEntries)
{
    neurons     = nullptr;
    synapses    = nullptr;
    recorder    = nullptr;
    stimulus    = nullptr;

    info.dt                         = 0.0;
    info.simulationTime             = 0.0;
	info.density					= 0;
	info.Lx							= 0.0;
	info.Ly							= 0.0;
	info.Dimensions					= 0;
    info.networkScaling_synStrength = 0.0;
    info.time_step                  = 0;
    info.networkScaling_extern      = 0;
    info.N                          = 0;
    info.globalSeed                 = -1;
    info.pathTo_inputFile 			= "";

    LoadParameters(baseDir,parEntries);

    SaveParameterOptions();
}



void NeuralNetwork::SaveParameters(){

    // if this test file does not appear in the target directory: stop the
    // simulation and check the directoryPath.

    std::ofstream stream(recorder->GetParametersFilename());
    recorder->WriteHeader(&stream);
    //stream <<  "#*****************************************************************\n";
    //stream <<  "#Date and Time:             " << std::ctime(&end_time);
    //stream <<  "#*****************************************************************\n";
    stream <<  "Title                       " << this->recorder->GetTitle()  << "\n";
    stream <<  "#*****************************************************************\n";
    stream <<  "simulationTime              " << std::to_string(info.simulationTime)  << " \t\tseconds\n";
    stream <<  "dt_timestep                 " << std::to_string(info.dt)  << " \t\tseconds\n";
	stream <<  "globalSeed                  " << std::to_string(info.globalSeed)  << " \t\t\t\t#overrides all other seeds if unequal -1\n";

	stream << "#*****************************************************************\n";
	stream << "#****************        Spatial parameters        ***************\n";
	stream << "#*****************************************************************\n";
	stream << "density                     " << std::to_string(info.density) << " \t\t\t\t#total number of neurons/mm^2 or /mm depending on Dimensions \n";
	stream << "Dimensions                  " << std::to_string(info.Dimensions) << " \t\t\t\t#for 1D set 1; for 2D set 2 \n";

    stream <<  "#*****************************************************************\n";
    stream <<  "#********** Scaling of synaptic and stimulus strengths ***********\n";
    stream <<  "#*****************************************************************\n";
    stream <<  "scalingSynapticStrength     " << std::to_string(info.networkScaling_synStrength)  << "\t\t#Scaling exponent. Set = 0 if no scaling needed, otherwise typical exponent is -0.5.\n";
    stream <<  "scaling_C_N                 " << std::to_string(info.networkScaling_extern)  << "\t\t\t\t# Set = 0 to scale with number of presynaptic neurons C. Set = 1 to scale with total number of neurons N. (details below) \n";
    stream <<  "#\t\tscaling_C_N=0    scales internal synaptic strengths and UncorrelatedStimulus with C^s    and    WhiteNoiseStimulus and SpatialGaussianStimulus with 1 \n";
    stream <<  "#\t\tscaling_C_N=1    scales internal synaptic strengths and UncorrelatedStimulus with N^s    and    WhiteNoiseRescaled and SpatialGaussianStimulus with N^(-s) \n";
    stream <<  "#\t\tscalingSynapticStrength = s, N = number of neurons from all populations, C = average number of presynaptic neurons.\n";


    this->neurons->SaveParameters(&stream);
    this->stimulus->SaveParameters(&stream);
    this->recorder->SaveParameters(&stream);
    this->synapses->SaveParameters(&stream);

    stream << std::endl;
    stream.close();
}

void NeuralNetwork::SaveParameterOptions(){// This function should have stuff moved to the heap (and deleted at the end of the function)
	//return;//SaveParameteroptions is disabled
	std::cout << recorder->GetParameterOptionsFilename()<<"\n";
    std::ofstream stream(recorder->GetParameterOptionsFilename());


	GlobalSimInfo mockInfo = info;
	mockInfo.globalSeed = 1;//because many of the classes do not allow a negative seed, the ParameterOptions file is generated with an adapted globalInfo


    stream <<  "#************************************************************************************************\n";
    stream <<  "#*************  Stimulus options     ************************************************************\n";
    stream <<  "#************************************************************************************************\n";

    std::vector<std::string>  stim_strs;
	std::string	str;

	str = "stimulus_meanCurrent ";
	for (unsigned long p = 0; p < neurons->GetTotalPopulations(); p++) {str += "5.0 ";}
	stim_strs.push_back(str + std::to_string(info.simulationTime));
	str = "stimulus_sigmaCurrent ";
	for (unsigned long p = 0; p < neurons->GetTotalPopulations(); p++) { str += "1.0 "; }
	stim_strs.push_back(str + std::to_string(info.simulationTime));
    WhiteNoiseStimulus        stimulus_wn(neurons,&stim_strs,&mockInfo);
    stimulus_wn.SaveParameters(&stream);


    WhiteNoiseRescaled        stimulus_wnr(neurons,&stim_strs,&mockInfo);
    stimulus_wnr.SaveParameters(&stream);


	str = "UncorrelatedStimulus_stimulus_J_X ";
	for (unsigned long p = 0; p < neurons->GetTotalPopulations(); p++) { str = str + std::to_string((double(p)+1)/10)+ " "; }
	stim_strs.push_back(str);
    stim_strs.push_back("UncorrelatedStimulus_stimulus_step "+ std::to_string(info.simulationTime)+ " 5.00000");
    UncorrelatedPoissonLikeStimulus stimulus_upls(neurons,&stim_strs,&mockInfo);
    stimulus_upls.SaveParameters(&stream);


	std::vector<std::string>  stim_strs_WNL;
	str = "stimulus_meanCurrent ";
	for (unsigned int p = 0; p < neurons->GetTotalPopulations(); p++) { str += "0.0 "; }
	for (unsigned int p = 0; p < neurons->GetTotalPopulations(); p++) { str += "5.0 "; }
	stim_strs_WNL.push_back(str + "0.0 "+std::to_string(info.simulationTime));
	str = "stimulus_sigmaCurrent ";
	for (unsigned int p = 0; p < neurons->GetTotalPopulations(); p++) { str += "0.0 "; }
	for (unsigned int p = 0; p < neurons->GetTotalPopulations(); p++) { str += "1.0 "; }
	stim_strs_WNL.push_back(str + "0.0 " + std::to_string(info.simulationTime));
	WhiteNoiseLinear	stimulus_wnl(neurons, &stim_strs_WNL, &mockInfo);
	stimulus_wnl.SaveParameters(&stream);


	std::vector<std::string>  stim_strs_SGS;
	stim_strs_SGS.push_back("stimulus_NumberOfGaussians           2");
	stim_strs_SGS.push_back("stimulus_X_position                  0.250000	 0.750000");
	stim_strs_SGS.push_back("stimulus_Y_position                  0.250000	 0.750000");
	str = "stimulus_maxCurrent_1 ";
	for (unsigned int p = 0; p < neurons->GetTotalPopulations(); p++) { str += "5.0 "; }
	stim_strs_SGS.push_back(str + std::to_string(info.simulationTime));
	str = "stimulus_sigmaCurrent_t_1 ";
	for (unsigned int p = 0; p < neurons->GetTotalPopulations(); p++) { str += "0.01 "; }
	stim_strs_SGS.push_back(str + std::to_string(info.simulationTime));
	stim_strs_SGS.push_back("stimulus_sigmaCurrent_x_1          1.0 "+ std::to_string(info.simulationTime));
	str = "stimulus_maxCurrent_2 ";
	for (unsigned int p = 0; p < neurons->GetTotalPopulations(); p++) { str += "5.0 "; }
	stim_strs_SGS.push_back(str + std::to_string(info.simulationTime));
	str = "stimulus_sigmaCurrent_t_2 ";
	for (unsigned int p = 0; p < neurons->GetTotalPopulations(); p++) { str += "0.01 "; }
	stim_strs_SGS.push_back(str + std::to_string(info.simulationTime));
	stim_strs_SGS.push_back("stimulus_sigmaCurrent_x_2          1.0 " + std::to_string(info.simulationTime));
	str = "stimulus_Background_Noise ";
	for (unsigned int p = 0; p < neurons->GetTotalPopulations(); p++) { str += "1.00 "; }
	stim_strs_SGS.push_back(str + std::to_string(info.simulationTime));
	SpatialGaussianStimulus	stimulus_sgs(neurons, &stim_strs_SGS, &mockInfo);
	stimulus_sgs.SaveParameters(&stream);


	stim_strs.push_back("stimulus_noExternalNeurons   1");
	stim_strs.push_back("stimulus_PoissonTableEntries 10");
	str = "ExternalConnectivity_lengthscale ";
	for (unsigned int p = 0; p < neurons->GetTotalPopulations(); p++) { str += "0.25 "; }
	stim_strs.push_back(str);
	str = "ExternalConnectivity_PeakProbability ";
	for (unsigned int p = 0; p < neurons->GetTotalPopulations(); p++) { str += "0.8 "; }
	stim_strs.push_back(str);
	str = "stimulus_J_X ";
	for (unsigned int p = 0; p < neurons->GetTotalPopulations(); p++) { str = str + std::to_string((double(p) + 1) / 10) + " "; }
	stim_strs.push_back(str);
	SpatialPoissonStimulus stimulus_sps(neurons, &stim_strs, &mockInfo);
	stimulus_sps.SaveParameters(&stream);

    stream <<  "#************************************************************************************************\n";
    stream <<  "#*************  Neuron options *****************************************************************\n";
    stream <<  "#************************************************************************************************\n";

    LIFNeuronPop neuronLIF(&mockInfo,0);
    neuronLIF.SaveParameters(&stream);

    EIFNeuronPop neuronEIF(&mockInfo,0);
    neuronEIF.SaveParameters(&stream);

    QIFNeuronPop neuronQIF(&mockInfo,0);
    neuronQIF.SaveParameters(&stream);

    PoissonNeuronPop neuronPoisson(&mockInfo,0);
    neuronPoisson.SaveParameters(&stream);

    DictatNeuronPop neuronInput(&mockInfo, 0);
    neuronInput.SaveParameters(&stream);

    stream <<  "#************************************************************************************************\n";
    stream <<  "#*************  Connectivity options ************************************************************\n";
    stream <<  "#************************************************************************************************\n";

    std::string syn;
    syn = "synapses_" + std::to_string(0) + "_" + std::to_string(0) + "_";
    CurrentSynapse synapse_current(&neuronLIF,&neuronLIF,&mockInfo);

    RandomConnectivity       geom_rc(&synapse_current, &mockInfo);
    geom_rc.SaveParameters(&stream, syn);

    stream << "#************************************************\n";
    AdjacencyMatrixConnectivity       geom_am(&synapse_current, &mockInfo);
    geom_am.SaveParameters(&stream, syn);

    stream <<  "#************************************************\n";
    BinaryRandomConnectivity geom_brc(&synapse_current,&mockInfo);
    geom_brc.SaveParameters(&stream,syn);
	stream << "#************************************************\n";
	DistanceConnectivity geom_dc(&synapse_current, &mockInfo);
	geom_dc.SaveParameters(&stream, syn);

    stream <<  "#************************************************************************************************\n";
    stream <<  "#*************  Synapse options *****************************************************************\n";
    stream <<  "#************************************************************************************************\n";

    synapse_current.SaveParameters(&stream,syn);

    stream <<  "#************************************************\n";
    MongilloSynapse synapse_ms(&neuronLIF,&neuronLIF,&mockInfo);
    synapse_ms.SaveParameters(&stream,syn);
    stream <<  "#************************************************\n";
    MongilloSynapseContinuous synapse_msc(&neuronLIF,&neuronLIF,&mockInfo);
    synapse_msc.SaveParameters(&stream,syn);
    stream <<  "#************************************************\n";
    PRGSynapseContinuous synapse_prgc(&neuronLIF,&neuronLIF,&mockInfo);
    synapse_prgc.SaveParameters(&stream,syn);
	stream << "#************************************************\n";
	PowerLawSynapse synapse_pl(&neuronLIF, &neuronLIF, &mockInfo);
	synapse_pl.SaveParameters(&stream, syn);
    stream <<  "#************************************************\n";
    ExponentialCurrentSynapse synapse_ecs(&neuronLIF,&neuronLIF,&mockInfo);
    synapse_ecs.SaveParameters(&stream,syn);
    stream <<  "#************************************************\n";

    stream <<  "#************************************************************************************************\n";
    stream <<  "#*************  Heterosynaptic plasticity alternative classes ***********************************\n";
    stream <<  "#************************************************************************************************\n";
    stream <<  "#************************************************\n";
    HeteroLIFNeuronPop hhneuronLIF(&mockInfo,0);
    hhneuronLIF.SaveParameters(&stream);
    stream <<  "#************************************************\n";

    HeteroCurrentSynapse hhsynapse_current(&hhneuronLIF,&hhneuronLIF,&mockInfo);
    hhsynapse_current.SaveParameters(&stream,syn);
    stream <<  "#************************************************\n";
    stream <<  "#************************************************************************************************\n";
    stream <<  "#*************  Heterosynaptic dendrites ********************************************************\n";
    stream <<  "#************************************************************************************************\n";
    stream <<  "#************************************************\n";
    MonoDendriteSTDPTazerart tazden(&mockInfo);
    tazden.SaveParameters(&stream,"neurons_0");

    stream <<  "#*************  Branched dendrites **************************************************************\n";
    BranchedResourceHeteroSTDP brhSTDP(&mockInfo);
    brhSTDP.SaveParameters(&stream,"neurons_0");
//Instead of LIF, heteroLIF, put morphology options from Tazerart and HCS and HCP

    stream << std::endl;
    stream.close();
}


int NeuralNetwork::LoadParameters(std::string baseDir,std::vector<FileEntry> *parEntries){

    std::vector<std::string>        values;
    std::vector<std::string>        full_strs,neur_strs,syn_strs,stim_strs,rec_strs;
    std::string                     str_line,name,rec_type,stim_type,neurons_type;
    std::string                     line,title = "";
    bool                            timestep_bool{false};

    for (auto & parEntry : *parEntries) {

        name   = parEntry.name;
        values = parEntry.values;

        if(name.find("Title") != std::string::npos){
            if(values.empty())
                title = "default";
            else
                title = values.at(0);
        }
        else if((name.find("dt_timestep") != std::string::npos || name.find("dt") != std::string::npos)&& !timestep_bool){
            info.dt = std::stod(values.at(0));
            timestep_bool=true;
        }
        else if(name.find("simulationTime") != std::string::npos){
            info.simulationTime = std::stod(values.at(0));
        }
        else if(name.find("globalSeed") != std::string::npos){
            info.globalSeed = static_cast<int> (std::stod(values.at(0)));
            info.globalGenerator  = std::mt19937(info.globalSeed);
        }
		else if (name.find("density") != std::string::npos) {
            info.density = static_cast<int>(std::stod(values.at(0)));
		}
		else if (name.find("Dimensions") != std::string::npos) {
            info.Dimensions = static_cast<int>(std::stod(values.at(0)));
		}
        else if((name.find("synapticScaling") != std::string::npos) ||
                (name.find("scalingSynapticStrength") != std::string::npos)){
            info.networkScaling_synStrength = std::stod(values.at(0));
        }
        else if(name.find("scaling_C_N") != std::string::npos){
            info.networkScaling_extern = std::stoi(values.at(0));
        }
        else if(name.find("neurons_type") != std::string::npos){
            neurons_type = values.at(0);//Is this still relevant ?
        }
        else if(name.find("recorder_type") != std::string::npos){
            rec_type = values.at(0);
        }
        else if(name.find("stimulus_type") != std::string::npos){
            stim_type = values.at(0);
        }
        else if (name.find("pathToInputFile") != std::string::npos) {
            if (values.size() == 0) {
                info.pathTo_inputFile = "";
            }
            else {
                info.pathTo_inputFile = values.at(0);
            }
        }
    }
    info.pathTo_inputFile = info.pathTo_inputFile + title + "_";
    
    /*FilterStringEntries(parEntries,"neurons_",  &neur_strs);
    FilterStringEntries(parEntries,"synapses_", &syn_strs);
    FilterStringEntries(parEntries,"stimulus_", &stim_strs);
    FilterStringEntries(parEntries,"recorder_", &rec_strs);*/

    //assemble vector of lines for classes that are not yet adapted to ParameterEntry structs
    for (auto & parEntry : *parEntries) {
        line.clear();
        line = parEntry.name + " ";
        for(auto & val : parEntry.values)
            line.append(val+" ");
        full_strs.push_back(line);
    }

    FilterStringVector(&full_strs,"neurons_",  &neur_strs);
    FilterStringVector(&full_strs,"synapses_", &syn_strs);
    FilterStringVector(&full_strs,"stimulus_", &stim_strs);
    FilterStringVector(&full_strs,"recorder_", &rec_strs);

    /*if(neurons_type == str_LIFNeuron)
        this->neurons   = new LIFNeuronPop(&info);
    else if(neurons_type == str_EIFNeuron)
        this->neurons   = new EIFNeuronPop(&info);
    else if(neurons_type == str_QIFNeuron)
        this->neurons   = new QIFNeuronPop(&info);
    else{
        std::cout << "Neuron not defined";
        return 0;
     }*/
    neurons   = new NeuronPopSample(&neur_strs,&info);
    synapses  = new SynapseSample(neurons,&syn_strs,&info);

    if(stim_type == str_uncorrelatedStimulus)
        this->stimulus  = new UncorrelatedPoissonLikeStimulus(neurons,&stim_strs,&info);
    else if(stim_type == str_whitenoiseStimulus)
        this->stimulus  = new WhiteNoiseStimulus(neurons,&stim_strs,&info);
    else if(stim_type == str_whitenoiseRescaled)
        this->stimulus  = new WhiteNoiseRescaled(neurons,&stim_strs,&info);
    else if (stim_type == str_whitenoiseLinear)
		this->stimulus = new WhiteNoiseLinear(neurons, &stim_strs, &info);
	else if (stim_type == str_spatialgaussianStimulus)
		this->stimulus = new SpatialGaussianStimulus(neurons, &stim_strs, &info);
	else if (stim_type == str_spatialpoissonStimulus)
		this->stimulus = new SpatialPoissonStimulus(neurons, &stim_strs, &info);
    else{
        std::cout << "Stimulus not defined";
        return 0;
    }

    if(rec_type == str_advancedRecorder)
        this->recorder  = new AdvancedRecorder(neurons,synapses,stimulus,baseDir,&rec_strs,title,&info);
    else{
        std::cout << "Recorder not defined";
        return 0;
    }

    return 1;
}

int NeuralNetwork::WellDefined(){
    if((neurons == nullptr) || (synapses == nullptr) || (recorder == nullptr) || (stimulus == nullptr)){
        std::cout << "***************************************************\n";
        std::cout << "Neurons, Synapses, Recorder or Stimulus not defined\n";
        std::cout << "***************************************************\n";
        return 0;
    }

    return 1;
}

int NeuralNetwork::Simulate()
{
    if(!WellDefined())
        return 1;

    //******************************
    // Declarations & Initialization
    //******************************
    double   r;
    double   t_comp = 0;
    unsigned int      P = neurons->GetTotalPopulations();
    std::uniform_real_distribution<double> uniformDistribution (0.0,1.0);
    long      simSteps      = static_cast<long>(info.simulationTime/info.dt);  // number of simulation time steps
    // int      global_D_max = this->synapses->GetMaxD();          // get maximum delay across all synapses: size of waiting matrix DEPRECATED

    std::vector<std::vector<double>> synaptic_dV;
    synaptic_dV.resize(P);
    for(unsigned int p = 0; p < P; p++)
		synaptic_dV[p].resize(neurons->GetNeuronsPop(p),0);

    // Initialize waiting matrix for synaptic delays (safe to delete?) DEPRECATED
	/*
    std::vector<std::vector<std::vector<double>>> waiting_matrix (global_D_max+1, std::vector<std::vector<double>>(P, std::vector<double>(0)));
    for(int d_index = 0; d_index < waiting_matrix.size(); d_index++){
        for(int p_index = 0; p_index < waiting_matrix[0].size(); p_index++){
            waiting_matrix[d_index][p_index].resize(neurons->GetNeuronsPop(p_index),0);
        }
    }*/

    this->recorder->SetFilenameDate();
    SaveParameters();

    //*****************************************************
    // --------------- START OF THE SIMULATION ------------
    //*****************************************************
    info.time_step = 0;
    info.waiting_index = 0;
    clock_t begin  = clock();
    clock_t intermediate_time;

    this->synapses->ConnectNeurons();
    this->recorder->WriteConnectivity();
    this->recorder->WriteDistributionD();
    this->recorder->WriteDistributionJ();
    this->recorder->WriteDataHeader();

	std::cout << "\n Pandas start simulation : " << this->recorder->GetTitle() << "\n";
	auto start = std::chrono::high_resolution_clock::now();
//Remember remember
    NeuralNetwork::outputHeteroEvents();
//The 5th of November
    while(info.time_step <= simSteps)
    {
        for (unsigned int p{0}; p < P; p++) {
            for (unsigned long i{ 0 }; i < neurons->GetNeuronsPop(p); i++)
				synaptic_dV[p][i] = 0.0;
		}
        this->synapses->advect(&synaptic_dV);
        this->stimulus->Update(&synaptic_dV);
        this->neurons->advect(&synaptic_dV);
		this->synapses->Reset();
        this->recorder->Record(&synaptic_dV);


        if(info.time_step%(static_cast<long>(simSteps*0.01)) == 1){
            intermediate_time = clock();
            r = (static_cast<double>(info.time_step) / static_cast<double>(simSteps));
            t_comp = static_cast<double>(intermediate_time - begin) / CLOCKS_PER_SEC;
            std::cout << static_cast<int>(r*100) << "%  -- Comp. time: " << static_cast<int>(t_comp) << "/" << static_cast<int>(t_comp/r) << " sec. -- " << " \n";
        }

        info.time_step++;
    }

    this->recorder->Record(&synaptic_dV);
    this->recorder->writeFinalDataFile(t_comp);
	std::cout << "\nPandas end simulation : " <<this->recorder->GetTitle() << "\n";
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    std::cout << "\nTotal simulation time(s): " << duration.count() << "\n";

    NeuralNetwork::outputHeteroEvents();
//We should also make sure this content is stored in an output file .txt. This function is a console output for debugging.

    //*****************************************************
    // --------------- END OF THE SIMULATION ------------
    //*****************************************************

    return 0;
}

void NeuralNetwork::makeInputCopy(const std::string& inputFile) {
    this->recorder->makeInputCopy(inputFile);
}

void NeuralNetwork::outputHeteroEvents(){
        unsigned long potentiationCount {0}, depressionCount {0}, inBetweeners {0}, stable {0};
        bool printOut{ false };
    for (unsigned long popId{ 0 }; popId < this->neurons->GetTotalPopulations(); popId++) {
        if (neurons->GetPop(popId)->HasHeterosynapticPlasticity()){
            HeteroNeuronPop* pop = dynamic_cast<HeteroNeuronPop*>(neurons->GetPop(popId)); 
            //This dynamic cast happens only twice at runtime.
            for (unsigned long nId{ 0 }; nId < pop->GetNoNeurons(); nId++) {
                unsigned long synCount{ pop->GetSynapseCount(nId) };
                for (unsigned long sId{ 0 }; sId < synCount; ++sId) {
                    double w{ pop->GetWeight(nId, sId) };
                    if (w > 1.8) {
                        potentiationCount++;
                    }else if (w < 0.2) {
                        depressionCount++;
                    }else if (w > 0.8 && w < 1.2) {
                        stable++;
                    }else {
                        inBetweeners++;
                    }
                }
            }
        }
    }

if (printOut){    //std::ofstream stream(this->recorder->GetDirectoryPath(),  );
    std::cout << "Potentiation Count: " << potentiationCount << std::endl;
    std::cout << "Depression Count: " << depressionCount << std::endl;
    std::cout << "Stable Count: " << stable << std::endl;
    std::cout << "InBetweener Count: " << inBetweeners << std::endl;
}
}