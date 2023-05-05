
#include "AdvancedRecorder.hpp"


AdvancedRecorder::AdvancedRecorder(NeuronPopSample *ns, SynapseSample *syn, Stimulus *sti, std::string baseDir, std::vector<std::string> *input, std::string str_t, GlobalSimInfo * info) :Recorder(ns, syn, sti, baseDir, input, str_t, info)
{
	trackSynapses = false;

	Heatmap = 0;

	unsigned int P = neurons->GetTotalPopulations();

	noRasterPlotNeurons.resize(P, 0);
	notrackNeuronPotentials.resize(P, 0);
	CurrentContributions.resize(P);
    noTrackHeteroSynapsePerTrackedNeuron.resize(P, 0);

    currentBin.potential.resize(P);
	currentBin.spiker_ratio.resize(P);
	currentBin.externalCurrent.resize(P);
	currentBin.synapticState.resize(P);
	currentBin.synapticCurrents.resize(P);
	currentBin.totalCurrentSquared_mean.resize(P);
	currentBin.no_recordedSynapses.resize(P);
	currentBin.totalCurrent_mean_N.resize(P);

	for (unsigned int i = 0; i < P; i++)
	{
		currentBin.no_recordedSynapses[i].resize(P);
		currentBin.synapticState[i].resize(P);
		currentBin.synapticCurrents[i].resize(P);
		currentBin.totalCurrent_mean_N.at(i).resize(neurons->GetNeuronsPop(i), 0);
	}

	LoadParameters(input);
	CurrentContrBin.resize((static_cast<size_t>(P) + 1)* static_cast<size_t>(CurrentContributions.sum()));
	if (Heatmap != 0) {
		Densimap.resize(P);
		currentBin.Heatmap.resize(P);
		for (unsigned int pop = 0;pop < P;pop++) {
            size_t heatmapDimensions{static_cast<size_t>(pow(Heatmap, info->Dimensions))};
            currentBin.Heatmap[pop].resize(heatmapDimensions);
			Densimap[pop].resize(heatmapDimensions);
			if (info->Dimensions == 2) {
				for (unsigned long cell = 0;cell < neurons->GetNeuronsPop(pop);cell++)
					Densimap[pop][static_cast<size_t>(Heatmap*(floor(neurons->GetY_Pos(pop, cell) * Heatmap / info->Ly)) + floor(neurons->GetX_Pos(pop, cell) * Heatmap / info->Lx))] += 1;
			}
			else if (info->Dimensions == 1) {
				for (unsigned long cell = 0;cell < neurons->GetNeuronsPop(pop);cell++)
					Densimap[pop][static_cast<size_t>(floor(neurons->GetX_Pos(pop, cell) * Heatmap / info->Lx))] += 1;
			}
		}
	}
	if (trackSynapses) {
		for (unsigned i = 0; i < P; i++) {
			for (unsigned j = 0; j < P; j++)
				currentBin.synapticState[i][j].resize(this->synapses->GetNumberOfDataColumns(i, j));
		}
	}
}

void AdvancedRecorder::SaveParameters(std::ofstream * stream){
    Recorder::SaveParameters(stream);

    double  dt = info->dt;
    *stream << "recorder_binSize\t\t\t\t\t" << std::to_string(GetAveragingSteps()*dt)  << " #seconds \t#Bin size over which data saved in main recording data file is average over\n";
    *stream << "recorder_noRasterPlotNeurons\t\t\t";
    for(unsigned i = 0; i < noRasterPlotNeurons.size();i++) {
		*stream << std::to_string(noRasterPlotNeurons[i]) << " \t";
	}
	*stream << std::to_string((static_cast<double>(raster_t_0))*info->dt) << "\t";
    *stream << "#Record spike times of x neurons for (i-th column is x for the i-th population). The i+1-th column sets the initial recording time. If negative, records all neurons of pop\n";

    *stream << "recorder_notrackNeuronProfiles\t\t";
    for(unsigned i = 0; i < notrackNeuronPotentials.size();i++)
        *stream << std::to_string(notrackNeuronPotentials[i]) << " \t";
    *stream << "\t\t\t#Record currents and potentials at all time steps of the first x_p neurons, p = population index. [column 1: track #neurons in pop1, column 2: track #neurons in pop2, .. ]\n";

	*stream << "recorder_CurrentContributions\t\t\t";

	for (unsigned i = 0; i < CurrentContributions.size();i++)
		*stream << std::to_string(CurrentContributions[i])<<"\t";
	*stream << std::to_string((static_cast<double>(current_t_0))*info->dt) << "\t";
	*stream << "#Record the sources of input current to x neurons. (i-th column is x for the i-th population). The i+1-th column sets the initial recording time\n";

    *stream <<  "recorder_trackSynapses\t\t\t\t" << std::to_string(trackSynapses)  << "\t\t\t\t\t#Set = 1 to track averaged data from synapes, Set = 0 to ignore.\n";
	*stream <<  "recorder_Heatmap\t\t\t\t\t" << std::to_string(Heatmap) << "\t\t\t\t\t#Number of bins used to represent each dimension of the spatial domain in the firing rates Heatmap\n";
 
    *stream <<  "recorder_notrackHeteroSynapseProfiles\t";
    for (unsigned i = 0; i < noTrackHeteroSynapsePerTrackedNeuron.size();i++)
        *stream << std::to_string(noTrackHeteroSynapsePerTrackedNeuron[i])<< "\t";
    *stream <<std::to_string(heteroRecordingPerSteps)<< "\t\t#Number of synapses per neuron in recorder_notrackNeuronProfiles to record (per pop) and last number is the timesteps for each record call\n";

    *stream <<  "recorder_parsing\t\t\t\t\t";
    if (parserEnabled){
        *stream << "ON";
    } else {
        *stream << "OFF";
    }
    *stream << "\t\t#Enabling parsing of rasterplot data into spiketimes. ON vs OFF.\n";
}

void AdvancedRecorder::LoadParameters(std::vector<std::string> *input){
    std::string              name,token;
    std::vector<std::string> values;

    for(std::vector<std::string>::iterator it = (*input).begin(); it != (*input).end(); ++it) {

        SplitString(&(*it),&name,&values);

		if (name.find("recorder_binSize") != std::string::npos) {
            SetAveraging_steps(std::stod(values.at(0)));
        } else if (name.find("recorder_noRasterPlotNeurons") != std::string::npos) {
            SetNoRasterplotNeurons(&values);
		} else if (name.find("recorder_trackSynapses") != std::string::npos) {
            trackSynapses = std::stoi(values.at(0));
		} else if ((name.find("recorder_notrackNeuronPotentials") != std::string::npos) ||
			(name.find("recorder_notrackNeuronProfiles") != std::string::npos)) {
			SetNoTrackNeuronPotentials(&values);
		} else if (name.find("recorder_Heatmap") != std::string::npos) {
            Heatmap = std::stoi(values.at(0));
		} else if (name.find("recorder_CurrentContributions") != std::string::npos) {
            SetNoCurrentContribution(&values);
        } else if (name.find("recorder_notrackHeteroSynapseProfiles") != std::string::npos) {
            SetNoTrackHeteroSynapseProfilesPerTrackedNeuronPerPop(&values);
		} else if (name.find("recorder_parsing") != std::string::npos){
            if (values.at(0).find("ON") != std::string::npos){
                parserEnabled=true;
            }
        }
    }

    Recorder::LoadParameters(input);
}

void AdvancedRecorder::SetNoTrackNeuronPotentials(std::vector<std::string> *values){
    unsigned int P = neurons->GetTotalPopulations();
    for(unsigned int i = 0; i < min_(P, values->size());i++){
        notrackNeuronPotentials[i] = std::stoul(values->at(i));
        if((notrackNeuronPotentials[i] >= neurons->GetNeuronsPop(i)) ||
            (notrackNeuronPotentials[i] < 0)){
            std::cout << "Potentials: Tracking all neurons of population " << i << "\n";
            notrackNeuronPotentials[i] = neurons->GetNeuronsPop(i);
        }
    }
}

void AdvancedRecorder::SetNoCurrentContribution(std::vector<std::string> *values) {
	unsigned int P = neurons->GetTotalPopulations();
	for (unsigned int i = 0;i < min_(P, values->size());i++) {
		CurrentContributions[i] = min_(neurons->GetNeuronsPop(i), std::stoul(values->at(i)));
		if (CurrentContributions[i] < 0) {
			std::cout << "Current contribution of all neurons in population " << i << "\n";
			CurrentContributions[i] = neurons->GetNeuronsPop(i);
		}
	}
	if (values->size() > P) {
		if (values->at(P)[0] != '#')// The comments in the code are not meant to be read upon execution
			current_t_0 = static_cast<int>(std::round(std::stod(values->at(P)) / info->dt));
	}
}

void AdvancedRecorder::SetNoRasterplotNeurons(std::vector<std::string> *values){
    unsigned int P = neurons->GetTotalPopulations();
    for(unsigned int i = 0; i < min_(P, values->size());i++){
        noRasterPlotNeurons[i] = std::stoul(values->at(i));
        if((noRasterPlotNeurons[i]  > neurons->GetNeuronsPop(i)) ||
           (noRasterPlotNeurons[i] < 0)){
            std::cout << "Rasterplot: Tracking all neurons of population "<< i << "\n";
            noRasterPlotNeurons[i]  = neurons->GetNeuronsPop(i);
        }
    }
	if (values->size() > P) {
		if (values->at(P)[0] != '#')
			raster_t_0 = static_cast<int>(std::round(std::stod(values->at(P)) / info->dt));
	}
}


void AdvancedRecorder::SetNoTrackHeteroSynapseProfilesPerTrackedNeuronPerPop(std::vector<std::string> *values) {
    int P = static_cast<int>(neurons->GetTotalPopulations());
    for(int i = 0; i < min_(P,static_cast<int>(values->size()));i++){//min() only makes sense if you remove the hash
        if((noTrackHeteroSynapsePerTrackedNeuron[i]  > neurons->GetPop(i)->GetNumberOfSynapses()) || (noTrackHeteroSynapsePerTrackedNeuron[i] < 0)){
                noTrackHeteroSynapsePerTrackedNeuron[i]=neurons->GetPop(i)->GetNumberOfSynapses();
        }
        noTrackHeteroSynapsePerTrackedNeuron[i] = std::stoi(values->at(i));
    }
    if (values->size() > P) {//This requires an extra number to record??
        try {
            this->heteroRecordingPerSteps = std::stoul(values->at(P));
        } catch (...) {
        }
    } //else heteroRecordingPerSteps is 1
}


void AdvancedRecorder::WriteDataHeader_Heatmap() {
	if (Heatmap == 0)
		return;
	std::ofstream file;
	int P = this->neurons->GetTotalPopulations();
	int dim = info->Dimensions;
	if (dim != 1 && dim != 2) {
		std::cout << "Heatmap can only be used in 1D or 2D\n";
		throw "Heatmap can only be used in 1D or 2D\n";
	}


	for (int pop = 0;pop < P;pop++) {
		file.open(this->GetHeatmapFilename() + std::to_string(pop)+".dat", std::ofstream::out | std::ofstream::trunc);
		WriteHeader(&file);
		file << "# Population:" + std::to_string(pop) + "\n";
		file << "# Dimension:" + std::to_string(dim) + "\n";
		file << "# L=" + std::to_string(info->Lx) + "\n";
		file << "#************************************\n";
		if (dim == 1) {
			file << "# column 1 : t (secs.) \n";
			file << "# columnn n+1: r_" + std::to_string(pop) + "(Hz) for neurons in position:\n";
			file << "#\t\t (n-1) * " + std::to_string(info->Lx / Heatmap) + "  <  x  <  n * " + std::to_string(info->Lx / Heatmap) + "\n";
		}
		else {
			file << "# column 1 : t (secs.) \n";
			file << "# columnn n+1: r_" + std::to_string(pop) + "(Hz) for neurons in position:\n";
			file << "#\t\t (n-1)%" + std::to_string(Heatmap) + " * " + std::to_string(info->Lx / Heatmap) + "  <  x  <  n%" + std::to_string(Heatmap) + " * " + std::to_string(info->Lx / Heatmap) + "\n";
			file << "#\t\t (n-1)//" + std::to_string(Heatmap) + " * " + std::to_string(info->Lx / Heatmap) + "  <  y  <  n//" + std::to_string(Heatmap) + " * " + std::to_string(info->Lx / Heatmap) + "\n";
		}
		//legend header
		if (dim == 1) {
			file << "t\t";
			for (int x = 0;x < Heatmap;x++) {
				file <<"["<< double(x) / double(Heatmap)*info->Lx<<"]\t";
			}
		}
		else if (dim == 2) {
			file << "t\t";
			for (int x = 0;x < Heatmap;x++) {
				for (int y = 0;y < Heatmap;y++) {
					file << "[" << double(x) / double(Heatmap)*info->Lx <<","<< double(y) / double(Heatmap)*info->Lx << "]\t";
				}
			}
		}
		file << "\n#************************************\n";
		for (int k = 1;k < 2 + pow(Heatmap, dim);k++)
			file << "#" + std::to_string(k) << "  \t";
		file << "\n";
		file.close();
	}
}

void AdvancedRecorder::WriteDataHeader_Averages(){

    int P = this->neurons->GetTotalPopulations();
    int col = 1;

    this->FileStreams.averagesFileStream.open (this->GetDataFilename(), std::ofstream::out | std::ofstream::trunc);
    WriteHeader(&this->FileStreams.averagesFileStream);
    this->FileStreams.averagesFileStream << "#************************************\n";
    this->FileStreams.averagesFileStream << "#Columns: \n";
    this->FileStreams.averagesFileStream << "#" +std::to_string(col)+ " t (secs.) \n"; col++;
    for(int i = 0; i < P; i++){
        this->FileStreams.averagesFileStream << "#" +std::to_string(col)+ " V_" << i << " (mV)\n";col++;}

    for(int i = 0; i < P; i++){
        this->FileStreams.averagesFileStream << "#" +std::to_string(col)+ " r_" << i << " (Hz)\n"; col++;}

    for(int i = 0; i < P; i++){
        this->FileStreams.averagesFileStream << "#" +std::to_string(col)+ " mu_Ext_" << i << "/tau_m (dmV/sec) \n"; col++;}

    for(int i = 0; i < P; i++)
    {
        for(int j = 0; j < P; j++){
            this->FileStreams.averagesFileStream << "#" +std::to_string(col)+ " mu_" << i << "_" << j << "/tau_m (dmV/sec) \n"; col++;}
    }
    for(int i = 0; i < P; i++){
        this->FileStreams.averagesFileStream << "#" +std::to_string(col)+ " mu_" << i << "/tau_m (dmV/sec) \n"; col++;}

    for(int i = 0; i < P; i++){
        this->FileStreams.averagesFileStream << "#" +std::to_string(col)+ " Quenched fluctuations of mu_" << i << "/tau_m (dmV/sec) = sqrt(PopAver(TempAver(mu_i)^2) - PopAver(TempAver(mu_i))^2) \n"; col++;}

    for(int i = 0; i < P; i++){
        this->FileStreams.averagesFileStream << "#" +std::to_string(col)+ " Temporal fluctuations of mu_" << i << "/tau_m*sqrt(dt) (dmV/sqrt(sec)) = sqrt(dt)*sqrt(PopAver(TempAver(mu_i^2) - TempAver(mu_i)^2))\n"; col++;}

	this->FileStreams.averagesFileStream << "t\t";
	for (int i = 0; i < P; i++)
		this->FileStreams.averagesFileStream << "V_" << i<<"\t";
	for (int i = 0; i < P; i++)
		this->FileStreams.averagesFileStream << "r_" << i << "\t";
	for (int i = 0; i < P; i++)
		this->FileStreams.averagesFileStream << "mu_ext_" << i << "\t";
	for (int i = 0; i < P; i++) {
		for (int j=0;j<P;j++)
			this->FileStreams.averagesFileStream << "mu_" << i << "_"<< j << "\t";
	}
	for (int i = 0; i < P; i++)
		this->FileStreams.averagesFileStream << "mu_tot_" << i << "\t";
	for (int i = 0; i < P; i++)
		this->FileStreams.averagesFileStream << "sigma_quenched_" << i << "\t";
	for (int i = 0; i < P; i++)
		this->FileStreams.averagesFileStream << "sigma_t_" << i << "\t";
	this->FileStreams.averagesFileStream << "\n#************************************\n";

    for(int k = 1;k<col;k++)
        this->FileStreams.averagesFileStream << "#"+std::to_string(k) << "\t\t";
    this->FileStreams.averagesFileStream << "\n";
}

void AdvancedRecorder::WriteDataHeader_Rasterplot(){

    if(noRasterPlotNeurons.sum() == 0)
        return;

    this->FileStreams.rasterplotFileStream.open (GetRasterplotFilename(), std::ofstream::out | std::ofstream::trunc);
    WriteHeader(&this->FileStreams.rasterplotFileStream);
    this->FileStreams.rasterplotFileStream << "#1 t (secs.) \t #2 neuron_id \t #3 neuron_pop_id \t\n";
	this->FileStreams.rasterplotFileStream << "# Note that the recorded neurons are equidistant within the population. Therefore, neuron_id are not necessarily successive numbers\n";
	this->FileStreams.rasterplotFileStream << "Spike_t\tNeuron_id\tPop_id\n";
	this->FileStreams.rasterplotFileStream << "#************************************\n";
}

void AdvancedRecorder::WriteDataHeader_Potential(){

    if(notrackNeuronPotentials.sum() == 0)
        return;

    unsigned long             P = neurons->GetTotalPopulations();
    this->FileStreams.potentialFileStream.open(GetPotentialFilename(), std::ofstream::out | std::ofstream::trunc);

    WriteHeader(&this->FileStreams.potentialFileStream);
    this->FileStreams.potentialFileStream << "#1 t (secs.)\t 2-"<<1+notrackNeuronPotentials.sum()<<" V_pop_id (mV) \n";

	this->FileStreams.potentialFileStream << "t\t";
    for(unsigned long p = 0;p<P;p++){
        for(unsigned long i = 0;i<notrackNeuronPotentials[p];i++)
            this->FileStreams.potentialFileStream << "V_" << p << "_" << i <<  "\t";
    }
	this->FileStreams.potentialFileStream << "\n#************************************\n";

}

void AdvancedRecorder::WriteDataHeader_Currents(){

    if(notrackNeuronPotentials.sum() == 0)
        return;

    unsigned long             P = neurons->GetTotalPopulations();
    this->FileStreams.currentsFileStream.open(GetCurrentsFilename(), std::ofstream::out | std::ofstream::trunc);

    WriteHeader(&this->FileStreams.currentsFileStream);
    this->FileStreams.currentsFileStream << "#1 t (sec)\t 2-"<<1+ notrackNeuronPotentials.sum()<<" mu_pop_id / tau_m(dmV / sec)\n";

	this->FileStreams.currentsFileStream << "t\t";
	for(unsigned long p = 0;p<P;p++){
        for(unsigned long i = 0;i<notrackNeuronPotentials[p];i++)
            this->FileStreams.currentsFileStream << "mu_" << p << "_" << i << "\t";
    }
	this->FileStreams.currentsFileStream << "\n#************************************\n";

}

void AdvancedRecorder::WriteDataHeader_CurrentsContribution() {
	if (CurrentContributions.sum() == 0)
		return;
	unsigned long P = neurons->GetTotalPopulations();
	long index;
	long HeaderIndex;
	long N;
	this->FileStreams.cCurrentsFileStream.open(GetCurrentCrontributionFilename(), std::ofstream::out | std::ofstream::trunc);

	WriteHeader(&this->FileStreams.cCurrentsFileStream);
	this->FileStreams.cCurrentsFileStream << "# 1 t (sec)\n";
	HeaderIndex = 2;
	for (unsigned long p_target = 0;p_target < P;p_target++) {
		N = CurrentContributions[p_target];
		if (N > 0) {
			for (unsigned long p_source = 0;p_source < P;p_source++) {
				this->FileStreams.cCurrentsFileStream <<"# "<< HeaderIndex <<"-"<< HeaderIndex +N -1 <<" Irecurrent_pop"<< p_target<<"<-pop"<< p_source<<"_neuronid / tau_m (mV/s)\n";
				HeaderIndex = HeaderIndex + N;
			}
			this->FileStreams.cCurrentsFileStream << "# "<<HeaderIndex << "-" << HeaderIndex + N -1 << " Ifeedforward_pop"<< p_target <<"_id / tau_m (mV/s)\n";
			HeaderIndex = HeaderIndex + N;
		}
	}
	this->FileStreams.cCurrentsFileStream << "# Note that the recorded neurons are equidistant within the population. Therefore, neuron_id are not necessarily successive numbers\n";

	this->FileStreams.cCurrentsFileStream << "t\t";
	for (unsigned long p_target = 0;p_target < P;p_target++) {
		for (unsigned long p_source = 0;p_source < P;p_source++) {
			for (long i = 0;i < CurrentContributions[p_target];i++) {
                //This is an integer type division, it does not need floor()
				index = i*neurons->GetNeuronsPop(p_target) / CurrentContributions[p_target];
				this->FileStreams.cCurrentsFileStream << "IRec" << p_target << "_" << p_source << "_" << index << "\t";
			}
		}
		for (long i = 0;i < CurrentContributions[p_target];i++) {
            //This is an integer type division, it does not need floor()
			index = i*neurons->GetNeuronsPop(p_target) / CurrentContributions[p_target];
			this->FileStreams.cCurrentsFileStream << "Iffd" << p_target << "_" << index << "\t";
		}
	}
	this->FileStreams.cCurrentsFileStream << "\n#************************************\n";
}


void AdvancedRecorder::WriteDataHeader_SynapseStates() {
    if(!trackSynapses)
        return;

    this->FileStreams.synStatesFileStream.open (this->GetSynapseStateFilename(), std::ofstream::out | std::ofstream::trunc);
    WriteHeader(&this->FileStreams.synStatesFileStream);
    this->FileStreams.synStatesFileStream << "#************************************\n";
    this->FileStreams.synStatesFileStream << "#Columns: \n";
    this->FileStreams.synStatesFileStream << "#1 t (secs.) \n";
	this->FileStreams.synStatesFileStream << synapses->GetDataHeader(2);
	this->FileStreams.synStatesFileStream << "t\t"<<synapses->GetUnhashedDataHeader() << "\n";
    this->FileStreams.synStatesFileStream << "# Attention, Synaptic state data is spike-induced : at each time step, only synapses from which the presynaptic neuron has spiked are measured. Data only tested for CurrentSynapse and MongilloSynapse\n";
    this->FileStreams.synStatesFileStream << "#************************************\n";
    return;

}


void AdvancedRecorder::WriteDataHeader_HeteroSynapses(){

    if(noTrackHeteroSynapsePerTrackedNeuron.sum() == 0)
        return;

    unsigned long P = neurons->GetTotalPopulations();
    unsigned long synTrackCount;

    this->FileStreams.heteroSynapsesFileStream.open(GetHeteroSynapseStateFilename(), std::ofstream::out | std::ofstream::trunc);
    WriteHeader(&this->FileStreams.heteroSynapsesFileStream);

    for(unsigned long p = 0;p<P;p++){
        synTrackCount =  noTrackHeteroSynapsePerTrackedNeuron[p];
        if (synTrackCount != 0 && this->neurons->GetPop(p)->HasHeterosynapticPlasticity()){
            this->FileStreams.heteroSynapsesFileStream << "#Pop. "<< p << " profile -> "<<this->neurons->GetPop(p)->GetIndividualSynapticProfileHeaderInfo() <<" \n";
        }
    }

    this->FileStreams.heteroSynapsesFileStream << "\n#************************************\n";
    this->FileStreams.heteroSynapsesFileStream << "#1 t (secs.)\t 2-"<<1+noTrackHeteroSynapsePerTrackedNeuron.sum()<<" Profile_pop_id_neuron_id_synapse_id \n";
    this->FileStreams.heteroSynapsesFileStream << "t\t";

    for(unsigned long p = 0;p<P;p++){
        synTrackCount =  noTrackHeteroSynapsePerTrackedNeuron[p];
        if ( synTrackCount == 0 || !this->neurons->GetPop(p)->HasHeterosynapticPlasticity()) {
            continue;
        }
        for(unsigned long i = 0;i<notrackNeuronPotentials[p];i++) {
            for (unsigned long k = 0; k < noTrackHeteroSynapsePerTrackedNeuron[p]; ++k) {
                this->FileStreams.heteroSynapsesFileStream << "Profile_" << p << "_" << i << "_" << (k) <<  "\t";
            }
        }
    }
    this->FileStreams.heteroSynapsesFileStream << "\n#************************************\n";
}

void AdvancedRecorder::WriteDataHeader_HeteroSynapsesOverall(){

    if(noTrackHeteroSynapsePerTrackedNeuron.sum() == 0)
        return;

    //std::cout << "The file has been properly created!!!!\n";
    unsigned long P = neurons->GetTotalPopulations();
    unsigned long synTrackCount;
    this->FileStreams.hSOverallFileStream.open(GetOverallHeteroSynapseStateFilename(), std::ofstream::out | std::ofstream::trunc);

    WriteHeader(&this->FileStreams.hSOverallFileStream);
    for(unsigned long p = 0;p<P;p++){
        synTrackCount =  noTrackHeteroSynapsePerTrackedNeuron[p];
        if (synTrackCount != 0 && this->neurons->GetPop(p)->HasHeterosynapticPlasticity()){
            this->FileStreams.heteroSynapsesFileStream << "#Pop. "<< p << " Overall Profile -> "<<this->neurons->GetPop(p)->GetOverallSynapticProfileHeaderInfo() <<" \n";
        }
    }
    this->FileStreams.hSOverallFileStream << "\n#************************************\n";

    this->FileStreams.hSOverallFileStream << "#1 t (secs.)\t 2-"<<1+noTrackHeteroSynapsePerTrackedNeuron.sum()<<" Profile_pop_id_neuron_id \n";

    this->FileStreams.hSOverallFileStream << "t\t";

    for(unsigned long p = 0;p<P;p++){
        synTrackCount =  noTrackHeteroSynapsePerTrackedNeuron[p];
        if (!this->neurons->GetPop(p)->HasHeterosynapticPlasticity() || synTrackCount == 0) {
            continue;
        }
        for(unsigned long i = 0;i<notrackNeuronPotentials[p];i++) {
                this->FileStreams.hSOverallFileStream << "OverallProfile_" << p << "_" << i <<  "\t";
        }
    }

    this->FileStreams.hSOverallFileStream << "\n#************************************\n";
}


void AdvancedRecorder::WriteDataHeader(){
    WriteDataHeader_Rasterplot();
    WriteDataHeader_Averages();
    WriteDataHeader_SynapseStates();
    WriteDataHeader_Potential();
    WriteDataHeader_Currents();
	WriteDataHeader_Heatmap();
	WriteDataHeader_CurrentsContribution();
	WriteDataHeader_HeteroSynapses();
	WriteDataHeader_HeteroSynapsesOverall();
    reset_statistics();
}

void AdvancedRecorder::reset_statistics()
{
    int P = this->neurons->GetTotalPopulations();

    currentBin.potential         = 0.0;
    currentBin.spiker_ratio      = 0.0;
    currentBin.externalCurrent   = 0.0;
    currentBin.totalCurrentSquared_mean  = 0.0;

    for(int i = 0; i < P; i++)
    {
        std::fill(currentBin.totalCurrent_mean_N.at(i).begin(), currentBin.totalCurrent_mean_N.at(i).end(), 0);
        for(int j = 0; j < P; j++)
        {
            currentBin.synapticCurrents[i][j]    = 0.0;
            currentBin.no_recordedSynapses[i][j] = 0;
            currentBin.synapticState[i][j]       = 0.0;
        }
		if (Heatmap != 0) {
			for (int pixel = 0;pixel < pow(Heatmap, info->Dimensions);pixel++)
				currentBin.Heatmap[i][pixel] = 0.0;
		}
    }
}

void AdvancedRecorder::Record_Rasterplot(){

	if (noRasterPlotNeurons.sum() == 0 || info->time_step < raster_t_0)
        return;
    double           dt = info->dt;
    std::vector<long> *spiker;
    // int				first_id;
	long			Nrec;//number of neurons to record
	long			Ntot;//number of neurons in population

    for(unsigned int pop = 0;pop<neurons->GetTotalPopulations();pop++){
		Nrec = noRasterPlotNeurons[pop];
		if (Nrec == 0)
			continue;
        spiker   = neurons->GetSpikers(pop);
		Ntot = neurons->GetNeuronsPop(pop);

        for(auto const &spiker_id : (*spiker) ){
            if (spiker_id == std::round(static_cast<double>(Ntot) / Nrec * std::round(static_cast<double>(Nrec * spiker_id / Ntot))))
                this->FileStreams.rasterplotFileStream << static_cast<double>(info->time_step)*dt <<"\t" << spiker_id << "\t" << pop << "\n";
        }
    }
}

void AdvancedRecorder::Record_Currents(std::vector<std::vector<double>> * synaptic_dV){

    if(notrackNeuronPotentials.sum() == 0)
        return;
    double           dt = info->dt;
    double           t = double(info->time_step)*dt;
    int              P = neurons->GetTotalPopulations();

    SaveDoubleFile(&this->FileStreams.currentsFileStream,t,5);
    for(int p = 0;p<P;p++){
        for(unsigned long i = 0;i<notrackNeuronPotentials[p];i++){
            SaveDoubleFile(&this->FileStreams.currentsFileStream,synaptic_dV->at(p).at(i)/info->dt,5);
        }
    }
    this->FileStreams.currentsFileStream << "\n";
}

void AdvancedRecorder::Record_Potential(){

    if(notrackNeuronPotentials.sum() == 0)
        return;

    double           dt = info->dt;
    double           t = double(info->time_step)*dt;
    long             P = neurons->GetTotalPopulations();

    SaveDoubleFile(&this->FileStreams.potentialFileStream,t,5);
    for(int p = 0;p<P;p++){
        for(unsigned long i = 0;i<notrackNeuronPotentials[p];i++)
            SaveDoubleFile(&this->FileStreams.potentialFileStream,neurons->GetPotential(p,i),5);
    }
    this->FileStreams.potentialFileStream << "\n";
}

void AdvancedRecorder::Record_CurrentContributions(std::vector<std::vector<double>> * synaptic_dV) {

	if (CurrentContributions.sum() == 0 || info->time_step < current_t_0)
		return;

	double dt = info->dt;
	double t = double(info->time_step)*dt;
	unsigned P = neurons->GetTotalPopulations();
	long index;
	long IndexMemory = 0;


	for (unsigned int p_target = 0;p_target < P;p_target++) {
		for (unsigned int p_source = 0;p_source < P;p_source++) {
			for (long i = 0;i < CurrentContributions[p_target];i++) {
				index = i*neurons->GetNeuronsPop(p_target) / CurrentContributions[p_target];
                CurrentContrBin[static_cast<size_t>(i) + static_cast<size_t>(IndexMemory)] += synapses->GetRecurrentInput(p_target, p_source, index);
			}
			IndexMemory += CurrentContributions[p_target];
		}
		for (long i = 0;i < CurrentContributions[p_target];i++) {
			index = i*neurons->GetNeuronsPop(p_target) / CurrentContributions[p_target];
			CurrentContrBin[static_cast<size_t>(i) + static_cast<size_t>(IndexMemory)] += stimulus->GetSignalArray(p_target, index);
		}
		IndexMemory += CurrentContributions[p_target];
	}

	if ((info->time_step) % this->averaging_steps == 0) {
		SaveDoubleFile(&this->FileStreams.cCurrentsFileStream, t, 5);
		for (unsigned int i = 0; i < CurrentContrBin.size(); i++) {
			SaveDoubleFile(&this->FileStreams.cCurrentsFileStream, CurrentContrBin[i], 5);
			CurrentContrBin[i] = 0;
		}
		this->FileStreams.cCurrentsFileStream << "\n";
	}
}

void AdvancedRecorder::Record_SynapseStates(){

    double dt {info->dt};
    double value;
    unsigned int P{neurons->GetTotalPopulations()};

    if(!trackSynapses)
        return;

    if((info->time_step+1)%this->averaging_steps == 0)
    {
        SaveDoubleFile(&this->FileStreams.synStatesFileStream,static_cast<double>(info->time_step)*dt,4);

        //Record synaptic strengths
        for(unsigned post_population = 0; post_population < P; post_population++)
        {
            for(unsigned pre_population = 0; pre_population < P; pre_population++)
            {
                if(currentBin.no_recordedSynapses[post_population][pre_population] > 0)
                {
                    for(unsigned data_length = 0; data_length < currentBin.synapticState[post_population][pre_population].size(); data_length++)
                    {
                        value = currentBin.synapticState[post_population][pre_population][data_length] /static_cast<double>(currentBin.no_recordedSynapses[post_population][pre_population]);
                        SaveDoubleFile(&this->FileStreams.synStatesFileStream,value,6);
                    }
                }
                else
                {
                    for(unsigned data_length = 0; data_length < currentBin.synapticState[post_population][pre_population].size(); data_length++)
                    {
                        value = NAN;
                        SaveDoubleFile(&this->FileStreams.synStatesFileStream,value,6);
                    }
                }
            }
        }
        this->FileStreams.synStatesFileStream << "\n";
    }
}

void AdvancedRecorder::Record_Averages(){

    double value;
    unsigned int    pops       = neurons->GetTotalPopulations();
    double dt          {info->dt};
    double deltaT_rec {static_cast<double>(this->averaging_steps)*dt};
    double n_aver      {static_cast<double>(averaging_steps)};
    double n;

    std::valarray<double> means(pops);
    std::valarray<double> popAver_SquaredMeanTime(pops);

    means        = 0.0; popAver_SquaredMeanTime = 0.0;

    for(unsigned int pop = 0; pop < pops; pop++){
        n = neurons->GetNeuronsPop(pop);
        for(unsigned long j = 0;j<neurons->GetNeuronsPop(pop);j++){
            means[pop]                   += currentBin.totalCurrent_mean_N.at(pop).at(j)/n/n_aver;
            popAver_SquaredMeanTime[pop] += pow(currentBin.totalCurrent_mean_N.at(pop).at(j)/n_aver,2)/n; // population average over square of temporal average [(I_alpha^i)^2]
        }
    }


    if((info->time_step)%this->averaging_steps == 0)
    {
        SaveDoubleFile(&this->FileStreams.averagesFileStream,static_cast<double>(info->time_step)*dt,4);

        //Record average potential
        for(unsigned int i = 0; i < pops; i++)
        {
            value      = currentBin.potential[i]/n_aver;
            SaveDoubleFile(&this->FileStreams.averagesFileStream,value,3);
        }
        //Record firing rate
        for(unsigned int i = 0; i < pops; i++)
        {
            value = currentBin.spiker_ratio[i]/deltaT_rec;
            SaveDoubleFile(&this->FileStreams.averagesFileStream,value,3);
        }
        //Record external input 'current' (in mV/sec)
        for(unsigned int i = 0; i < pops; i++)
        {
            value   = currentBin.externalCurrent[i]/n_aver;
            SaveDoubleFile(&this->FileStreams.averagesFileStream,value,3);
        }
        //Record population averaged synaptic input currents mu_ij (in mV/sec)
        for(unsigned int i = 0; i < pops; i++)
        {
            for(unsigned int j = 0; j < pops; j++)
            {
                value   = currentBin.synapticCurrents[i][j]/n_aver;
                SaveDoubleFile(&this->FileStreams.averagesFileStream, value,3);
            }
        }
        //Record population averaged total synaptic input currents mu_i (in mV/sec)
        //for(int i = 0; i < pops; i++)
        //{
        //    value   = currentBin.totalCurrent_mean[i]/n_aver;
        //    SaveDoubleFile(&this->FileStreams.averagesFileStream,value,3);
        //}

        //Record check of mean total input currents
        for(unsigned int i = 0; i < pops; i++)
        {
            value = means[i];
            SaveDoubleFile(&this->FileStreams.averagesFileStream,value,3);
        }

        //Quenched fluctuations (due to varying input over population)
        for(unsigned int i = 0; i < pops; i++)
        {
            value = sqrt(popAver_SquaredMeanTime[i]-pow(means[i],2));
            SaveDoubleFile(&this->FileStreams.averagesFileStream,value,3);
        }

        //Temporal fluctuations
        for(unsigned int i = 0; i < pops; i++)
        {
            value = sqrt(currentBin.totalCurrentSquared_mean[i]/n_aver - popAver_SquaredMeanTime[i])*sqrt(info->dt);
            SaveDoubleFile(&this->FileStreams.averagesFileStream,value,3);
        }

        this->FileStreams.averagesFileStream << "\n";
        reset_statistics();
    }
}
//
void AdvancedRecorder::Record_Heatmap() {

	if (Heatmap == 0)
		return;

	std::ofstream file;
	double value;
	unsigned int    pops = neurons->GetTotalPopulations();
	double dt = info->dt;
	double deltaT_rec = double(this->averaging_steps)*dt;
	// double n_aver = double(averaging_steps);
	// double n;

	if ((info->time_step) % this->averaging_steps == 0){
		for (unsigned int pop = 0; pop < pops;pop++) {
			file.open(this->GetHeatmapFilename() + std::to_string(pop)+".dat", std::ofstream::out | std::ofstream::app);
			SaveDoubleFile(&file, double(info->time_step)*dt, 4);
			for (int pixel = 0;pixel < pow(Heatmap, info->Dimensions); pixel++) {
				value = currentBin.Heatmap[pop][pixel] / deltaT_rec;
				SaveDoubleFile(&file, value, 3);
			}
			file << "\n";
			file.close();
		}
	}
}

void AdvancedRecorder::Record_HeteroSynapses() {

    if(noTrackHeteroSynapsePerTrackedNeuron.sum() == 0) {
        return;
    }

    double           dt = info->dt;
    double           t = double(info->time_step)*dt;
    unsigned long             P = neurons->GetTotalPopulations();

    SaveDoubleFile(&this->FileStreams.heteroSynapsesFileStream,t,5);


    unsigned long synTrackCount;

    for(unsigned long p = 0;p<P;p++){
        synTrackCount =  noTrackHeteroSynapsePerTrackedNeuron[p];
        if (!this->neurons->GetPop(p)->HasHeterosynapticPlasticity() || synTrackCount == 0) {
            continue;
        }
        for(unsigned long i = 0;i<notrackNeuronPotentials[p];i++) {
            for (unsigned long k = 0; k < noTrackHeteroSynapsePerTrackedNeuron[p]; ++k) {
                SaveTupleOfDoublesFile(&this->FileStreams.heteroSynapsesFileStream, this->neurons->GetPop(p)->GetIndividualSynapticProfile(i, k), 5);
            }
        }
    }
    this->FileStreams.heteroSynapsesFileStream << "\n";
}

void AdvancedRecorder::Record_HeteroSynapsesOverall() {
    if(noTrackHeteroSynapsePerTrackedNeuron.sum() == 0) {
        return;
    }

    double           dt = info->dt;
    double           t = static_cast<double>(info->time_step)*dt;
    unsigned long             P = neurons->GetTotalPopulations();

    SaveDoubleFile(&this->FileStreams.hSOverallFileStream,t,5);

    unsigned long synTrackCount;

    for(unsigned long p = 0;p<P;p++){
        synTrackCount =  noTrackHeteroSynapsePerTrackedNeuron[p];
        if (!this->neurons->GetPop(p)->HasHeterosynapticPlasticity() || synTrackCount == 0) {
            continue;
        }
        for(unsigned long i = 0;i<notrackNeuronPotentials[p];i++) {
                SaveTupleOfDoublesFile(&this->FileStreams.hSOverallFileStream, this->neurons->GetPop(p)->GetOverallSynapticProfile(i), 5);
            //Here is selecting only the average weight per neuron, with precision 5 digits.
        }
    }
    this->FileStreams.hSOverallFileStream << "\n";
}

void AdvancedRecorder::Record(std::vector<std::vector<double>> * synaptic_dV)
{
    this->stepCount++;

	std::vector<long> *spikers;
	unsigned int P = this->neurons->GetTotalPopulations();
	double n;
	int Xbin, Ybin;

	for (unsigned int p = 0; p < P; p++) {
		for (unsigned long i = 0; i < neurons->GetNeuronsPop(p); i++) {
			n = static_cast<double>(this->neurons->GetNeuronsPop(p));
			currentBin.potential[p] += this->neurons->GetPotential(p, i) / n;
			currentBin.externalCurrent[p] += this->stimulus->GetSignalArray(p, i) / n / info->dt;
            //How is this arithmetic happening? Sequential division, first divide by n, then by dt.
			currentBin.totalCurrentSquared_mean[p] += pow(synaptic_dV->at(p).at(i) / info->dt, 2.0) / n;
			//record per neuron
			currentBin.totalCurrent_mean_N.at(p).at(i) += synaptic_dV->at(p).at(i) / info->dt;
		}
	}

	for (unsigned int prePop = 0;prePop < P; prePop++) {
		spikers = this->neurons->GetSpikers(prePop);
		currentBin.spiker_ratio[prePop] += static_cast<double>(spikers->size()) / static_cast<double>(this->neurons->GetNeuronsPop(prePop));
		if (Heatmap != 0) {
			if (info->Dimensions == 2) {
				for (unsigned int spike = 0;spike < spikers->size();spike++) {
					Xbin = static_cast<int>(neurons->GetX_Pos(prePop, spikers->at(spike)) * Heatmap / info->Lx);
					Ybin = static_cast<int>(neurons->GetY_Pos(prePop, spikers->at(spike)) * Heatmap / info->Ly);
                    currentBin.Heatmap[prePop][static_cast<size_t>(Heatmap) * Ybin + Xbin] += 1 / static_cast<double>(Densimap[prePop][static_cast<size_t>(Heatmap) * Ybin + Xbin]);
				}
			}
			else {
				for (unsigned int spike = 0;spike < spikers->size();spike++) {
					Xbin = static_cast<int>(neurons->GetX_Pos(prePop, spikers->at(spike)) * Heatmap / info->Lx);
					currentBin.Heatmap[prePop][Xbin] += 1 / static_cast<double> (Densimap[prePop][Xbin]);
				}
			}
		}
		for (unsigned int postPop = 0; postPop < P; postPop++)
			currentBin.synapticCurrents[postPop][prePop] += this->synapses->GetCumulatedDV(postPop, prePop) / static_cast<double>(this->neurons->GetNeuronsPop(postPop)) / info->dt;

		if (!trackSynapses)
			continue;
		for (unsigned int i = 0; i < spikers->size(); i++)
		{
			for (unsigned int post_population = 0; post_population < P; post_population++) {
                //std::cout << "Synapse " << std::to_string(prePop) << " to " << std::to_string(post_population) << " at time " << std::to_string(this->info->time_step) << "\n";
				currentBin.synapticState[post_population][prePop] += this->synapses->GetSynapticState(post_population, prePop, (*spikers)[i]);  //arrays are added
				currentBin.no_recordedSynapses[post_population][prePop] += this->synapses->GetNumberOfPostsynapticTargets(post_population, prePop, (*spikers)[i]);
                for(unsigned int i_data = 0; i_data<currentBin.synapticState[post_population][prePop].size(); i_data++){
                //std::cout << std::to_string(currentBin.synapticState[post_population][prePop][i_data]) << " ";
                //std::cout << std::to_string(currentBin.synapticState[post_population][prePop][i_data]/currentBin.no_recordedSynapses[post_population][prePop]) << "\n";
                }

            }
		}
	}
	Record_SynapseStates();
	Record_Heatmap();
	Record_Rasterplot();
	Record_Currents(synaptic_dV);
	//Record_Histogram(synaptic_dV);
	Record_Averages();
	Record_Potential();
	Record_CurrentContributions(synaptic_dV);
    if (this->heteroRecordingPerSteps !=0 && (this->stepCount % this->heteroRecordingPerSteps) == 0) {
        //This condition will trigger if there is a "extra" number in the Parameters.txt and the step count is divisible by that number
        //1 for every step recording, 2 for every two steps, etc...
        Record_HeteroSynapses();
    	Record_HeteroSynapsesOverall();
    }
}

void AdvancedRecorder::CloseStreams()
{
    FileStreams.averagesFileStream.close();
    if (Heatmap!=0){
        FileStreams.heatmapFileStream.close();
    } if (noRasterPlotNeurons.sum()!=0){
        FileStreams.rasterplotFileStream.close();
    } if (notrackNeuronPotentials.sum()!=0){
        FileStreams.potentialFileStream.close();
        FileStreams.currentsFileStream.close();
    } if (CurrentContributions.sum()!=0){
        FileStreams.cCurrentsFileStream.close();
    } if (trackSynapses){
        FileStreams.synStatesFileStream.close();
    } if (noTrackHeteroSynapsePerTrackedNeuron.sum()!=0){
        FileStreams.heteroSynapsesFileStream.close();
        FileStreams.hSOverallFileStream.close();
    // } if (streamingNOutputBool){
    //     for (std::ofstream& stream : FileStreams.neuronOuputFileStreams){
    //         stream.close();
    // }
    }
}

void AdvancedRecorder::writeFinalDataFile(double comp_time)
{
    std::chrono::time_point<std::chrono::system_clock> t;
    t                    = std::chrono::system_clock::now();
    std::time_t end_time {std::chrono::system_clock::to_time_t(t)};

    std::ofstream stream(GetParametersFilename(),std::ofstream::out | std::ofstream::app);
    stream <<  "#*****************************************************************\n";
    stream <<  "#Comp. finalized: " << std::ctime(&end_time);
    stream <<  "#Comp. time: " << comp_time << " secs.\n";
    stream.close();
    std::cout << "Results saved in " << directoryPath << "\n";

}
