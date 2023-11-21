
#include "Recorder.hpp"
// extern char const *const GIT_COMMIT;
Recorder::Recorder(const std::shared_ptr<NeuronPopSample> neurons, const std::shared_ptr<SynapseSample> synapses, const std::shared_ptr<Stimulus> stimulus, std::string baseDirectory, std::vector<FileEntry> inputParameters, std::string titleString, std::string nonIterateTitle, GlobalSimInfo* infoGlobal)
 : infoGlobal{infoGlobal}, simulationTitle{titleString}, nonIterateTitle{nonIterateTitle},directoryPath{baseDirectory}, neurons{neurons}, synapses{synapses}, stimulus{stimulus} {

    PopInt totalNeuronPops = neurons->GetTotalPopulations();

	noRasterPlotNeurons.resize(totalNeuronPops, 0);
	neuronPotentialsToRecord.resize(totalNeuronPops, 0);
	currentContributionsToRecord.resize(totalNeuronPops);
    heteroSynTracker.resize(totalNeuronPops, std::pair<NeuronInt, signed long>(0,0));

    currentBin.potential.resize(totalNeuronPops, 0);
	currentBin.spikerRatio.resize(totalNeuronPops, 0);
	currentBin.externalCurrent.resize(totalNeuronPops, 0);
	currentBin.synapticState.resize(totalNeuronPops);
	currentBin.synapticCurrents.resize(totalNeuronPops);
	currentBin.totalCurrentSquared_mean.resize(totalNeuronPops, 0);
	currentBin.noRecordedSynapses.resize(totalNeuronPops);
	currentBin.neuronTotalCurrentMean.resize(totalNeuronPops);

	for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
		currentBin.noRecordedSynapses.at(neuronPop).resize(totalNeuronPops);
		currentBin.synapticState.at(neuronPop).resize(totalNeuronPops);
		currentBin.synapticCurrents.at(neuronPop).resize(totalNeuronPops);
		currentBin.neuronTotalCurrentMean.at(neuronPop).resize(neurons->GetNeuronsPop(neuronPop));
	}

	LoadParameters(inputParameters);
	currentContrBin.resize((static_cast<size_t>(totalNeuronPops) + 1)* static_cast<size_t>(std::reduce(currentContributionsToRecord.begin(), currentContributionsToRecord.end())));
	if (recordedHeatmap != 0) {
		densimap.resize(totalNeuronPops);
		currentBin.heatmap.resize(totalNeuronPops);
		for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
            size_t heatmapDimensions{static_cast<size_t>(pow(recordedHeatmap, infoGlobal->dimensions))};
            currentBin.heatmap.at(neuronPop).resize(heatmapDimensions, 0.0);
			densimap.at(neuronPop).resize(heatmapDimensions, 0.0);
			if (infoGlobal->dimensions == 2) {
				for (NeuronInt neuron : std::ranges::views::iota(0,neurons->GetNeuronsPop(neuronPop))){
					densimap.at(neuronPop).at(static_cast<size_t>(recordedHeatmap*(floor(neurons->GetYPosition(neuronPop, neuron) * recordedHeatmap / infoGlobal->yAxisLength)) + floor(neurons->GetXPosition(neuronPop, neuron) * recordedHeatmap / infoGlobal->xAxisLength))) += 1;
                }
			}
			else if (infoGlobal->dimensions == 1) {
				for (NeuronInt neuron : std::ranges::views::iota(0,neurons->GetNeuronsPop(neuronPop))){
					densimap.at(neuronPop).at(static_cast<size_t>(floor(neurons->GetXPosition(neuronPop, neuron) * recordedHeatmap / infoGlobal->xAxisLength))) += 1;
                }
			}
		}
	}
	if (trackSynapses) {
		for (PopInt targetPop : std::ranges::views::iota(0, totalNeuronPops)) {
			for (PopInt sourcePop : std::ranges::views::iota(0, totalNeuronPops)){
                if (synapses->GetConnectedState(targetPop, sourcePop)){
       				currentBin.synapticState.at(targetPop).at(sourcePop).resize(this->synapses->GetNoDataColumns(targetPop, sourcePop));
                }
            }
		}
	}
}


void Recorder::SetFilenameDate(){
	bool Windows = false;

    #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
        Windows = true;
    #endif
	// auto localTime = std::chrono::floor<std::chrono::seconds>(std::chrono::zoned_time{std::chrono::current_zone(), std::chrono::system_clock::now()}.get_local_time()); //Does not work in GCC
    time_t timeInType = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    // dateTime =  std::format("{:%Y_%m_%d_%H-%M-%S}", localTime); //Does not work in GCC
    // std::to_string(yearMonthYear.year()) + "_" +
    // std::to_string(yearMonthYear.year()) + "_" +
    // std::to_string(timePtr->tm_mday) + "_" +
    // std::to_string(timePtr->tm_hour) + "-" +
    // std::to_string(timePtr->tm_min) + "-" +
    // std::to_string(timePtr->tm_sec);
    std::stringstream outputString;
    // tm timeStruct = threadsafe::localtime(&timeInType);
    threadsafe::put_time(timeInType, "%Y_%m_%d_%H-%M-%S",outputString);
    dateTime = outputString.str();
	
    if (Windows)
        directoryPath += simulationTitle + "_" + dateTime + "\\";
    else
        directoryPath += simulationTitle + "_" + dateTime +  "/";

#ifdef _WIN32
    _mkdir(directoryPath.c_str()); //0744
#elif __APPLE__
    mkdir(directoryPath.c_str(),0744);
#elif __linux__
    mkdir(directoryPath.c_str(),0744);
#endif


}

void Recorder::WriteConnectivity() const{
    this->synapses->WriteConnectivity(GetConnectivityFilename(),noNeuronsConnectivity);
}

void Recorder::WriteDistributionD() const{
    this->synapses->WriteDistributionD(GetDelayFilename(),noNeuronsDelay);
}

void Recorder::WriteDistributionJ() const {
    this->synapses->WriteDistributionJ(GetJPotFilename(),noNeuronsJPot);
}

void Recorder::SaveParameters(std::ofstream& wParameterStream) const {
    wParameterStream <<  "#*************************************************\n";
    wParameterStream <<  "#************** Recorder Parameters **************\n";
    wParameterStream <<  "#*************************************************\n";
    //wParameterStream <<  "recorder_type\t\t\t\t\t\t" << GetType() <<  "\n";
    wParameterStream <<  "recorder_noNeuronsConnectivity\t\t" << std::to_string(noNeuronsConnectivity) << "\t\t\t#saves connectivity matrices for the first x neurons of each populations\n";
    wParameterStream <<  "recorder_noNeuronsDelay\t\t\t" << std::to_string(noNeuronsDelay) << "\t\t\t#saves delay connectivity matrices for the first x neurons of each populations\n";
    wParameterStream <<  "recorder_noNeuronsJPot\t\t\t" << std::to_string(noNeuronsJPot) << "\t\t\t\t#saves Jpot connectivity matrices for the first x neurons of each populations\n";
    //Adv Recorder SP
    wParameterStream << "recorder_binSize\t\t\t" << std::to_string(GetAveragingSteps()*infoGlobal->dtTimestep)  << " #secs\t\t#Bin size over which data saved in main recording data file is average over\n";
    wParameterStream << "recorder_noRasterPlotNeurons\t\t";
    for(PopInt neuronPop : std::ranges::views::iota (0, static_cast<PopInt>(noRasterPlotNeurons.size()))) {
		wParameterStream << std::to_string(noRasterPlotNeurons.at(neuronPop)) << " ";
	}
    wParameterStream << "\t"<< std::to_string((static_cast<double>(startRecordingTime))*infoGlobal->dtTimestep) ;
    wParameterStream << "\t\t#Record spike times of x neurons for (i-th column is x for the i-th population). The i+1-th column sets the initial recording time. If negative, records all neurons of neuronPop\n";

    wParameterStream << "recorder_notrackNeuronProfiles\t\t";
    for(PopInt neuronPop : std::ranges::views::iota(0, static_cast<PopInt>(neuronPotentialsToRecord.size()))){
        wParameterStream << std::to_string(neuronPotentialsToRecord.at(neuronPop)) << " ";
    }
    wParameterStream << "\t\t\t#Record currents and potentials at all time steps of the first x_p neurons, totalNeuronPops = population index. [column 1: track #neurons in pop1, column 2: track #neurons in pop2, .. ]\n";

	wParameterStream << "recorder_CurrentContributions\t\t";

	for (PopInt neuronPop : std::ranges::views::iota(0, static_cast<PopInt>(currentContributionsToRecord.size()))){
		wParameterStream << std::to_string(currentContributionsToRecord.at(neuronPop))<<" ";
    }
	wParameterStream << "\t"<< std::to_string((static_cast<double>(initialCurrent))*infoGlobal->dtTimestep) ;
	wParameterStream << "\t#Record the sources of input current to x neurons. (i-th column is x for the i-th population). The i+1-th column sets the initial recording time\n";

    wParameterStream <<  "recorder_trackSynapses\t\t\t" << std::to_string(trackSynapses)  << "\t\t\t#Set = 1 to track averaged data from synapes, Set = 0 to ignore.\n";
	wParameterStream <<  "recorder_Heatmap\t\t\t" << std::to_string(recordedHeatmap) << "\t\t\t#Number of bins used to represent each dimension of the spatial domain in the firing rates Heatmap\n";
 
    wParameterStream <<  "recorder_notrackHeteroSynapticProfiles\t";
    for (PopInt neuronPop : std::ranges::views::iota(0, static_cast<PopInt>(heteroSynTracker.size()))){
        wParameterStream << std::to_string(heteroSynTracker.at(neuronPop).first)<< " "<< std::to_string(heteroSynTracker.at(neuronPop).second)<< "  ";
    }
    wParameterStream <<std::to_string(heteroRecordPerSteps)<< "\t#Col1: number of neurons to track in pop 0, col2: number of synapses to track in pop0, ... Col2P: number of synapses to track in popP, Col2P+1: record every N steps (default 10)\n";

    wParameterStream <<  "recorder_parsing\t\t\t";
    if (parserEnabled){
        wParameterStream << "ON";
    } else {
        wParameterStream << "OFF";
    }
    wParameterStream << "\t\t\t#Enabling parsing of rasterplot data into spiketimes. ON vs OFF.\n";
}


void Recorder::LoadParameters(const std::vector<FileEntry>& recorderParameters){

    for(auto& [parameterName, parameterValues]: recorderParameters) {

        if(parameterName.find("recorder_noNeuronsConnectivity") != std::string::npos){
            noNeuronsConnectivity = std::stoi(parameterValues.at(0));
        } else if(parameterName.find("recorder_noNeuronsDelay") != std::string::npos){
            noNeuronsDelay = std::stol(parameterValues.at(0));
        } else if(parameterName.find("recorder_noNeuronsJPot") != std::string::npos){
            noNeuronsJPot = std::stol(parameterValues.at(0));
        //AdvRecorder LP
        } else if (parameterName.find("recorder_binSize") != std::string::npos) {
            SetAveragingSteps(std::stod(parameterValues.at(0)));
        } else if (parameterName.find("recorder_noRasterPlotNeurons") != std::string::npos) {
            SetNoRasterplotNeurons(parameterValues);
		} else if (parameterName.find("recorder_trackSynapses") != std::string::npos) {
            trackSynapses = std::stoi(parameterValues.at(0));
		} else if ((parameterName.find("recorder_notrackNeuronPotentials") != std::string::npos) ||
			(parameterName.find("recorder_notrackNeuronProfiles") != std::string::npos)) {
			SetNoRecordedNeuronPotentials(parameterValues);
		} else if (parameterName.find("recorder_Heatmap") != std::string::npos) {
            recordedHeatmap = std::stoi(parameterValues.at(0));
		} else if (parameterName.find("recorder_CurrentContributions") != std::string::npos) {
            SetNoCurrentContribution(parameterValues);
        } else if (parameterName.find("recorder_notrackHeteroSynapticProfiles") != std::string::npos) {
            SetNoRecordedHeteroSynapticProfilesPerTrackedNeuronPerPop(parameterValues);
		} else if (parameterName.find("recorder_parsing") != std::string::npos){
            if (parameterValues.at(0).find("ON") != std::string::npos){
                parserEnabled=true;
            }
        }
    }
    if(std::reduce(noRasterPlotNeurons.begin(), noRasterPlotNeurons.end()) != 0){
        recordRasterPlot=true;
    }
    if(std::reduce(neuronPotentialsToRecord.begin(), neuronPotentialsToRecord.end()) != 0){
        recordNeuronPotentials=true;
    }
    if(std::reduce(currentContributionsToRecord.begin(), currentContributionsToRecord.end()) != 0){
		recordCurrentContributions=true;
    }
    if(std::accumulate(heteroSynTracker.begin(), heteroSynTracker.end(), 0, [](int accumulator,std::pair<NeuronInt, signed long> pop){return accumulator + (pop.first*pop.second);}) != 0){
        recordHeteroSynapses=true;
    }
}


void Recorder::SetAveragingSteps(double secondsPerBin) {
    double dtTimestep = infoGlobal->dtTimestep;
    if(secondsPerBin < dtTimestep){
        timeStepsPerBin = 1;
    } else {
        timeStepsPerBin = static_cast<int>(std::round(secondsPerBin / dtTimestep));
    }
}

void Recorder::BindNoHeteroSynapsesPerPop(PopInt neuronPop) {
    if((heteroSynTracker.at(neuronPop).second  > neurons->GetPop(neuronPop)->GetNoSynapses()) || (heteroSynTracker.at(neuronPop).second < 0)){
            heteroSynTracker.at(neuronPop).second=neurons->GetPop(neuronPop)->GetNoSynapses();
    }
}

void Recorder::AllocateAndAssignStreamBuffer(std::ofstream &outputStream) {
    bufferPtrs.push_back(std::make_unique<std::vector<char>>(customBufferSize));
    outputStream.rdbuf()->pubsetbuf(bufferPtrs.back()->data(), customBufferSize);//Second argument is supposed to be the size of buffer in number of chars, NOT BITS
}

void Recorder::WriteHeader(std::ofstream& fileStream) const {
    fileStream <<  "#*****************************************************************\n";
    fileStream <<  "# Time and Title: " << dateTime << " -- " << simulationTitle << "\n";
    fileStream <<  "#*****************************************************************\n";
}

void Recorder::MakeInputCopies(const std::string& inputFileAddress) {
    //Parameters file
    std::ifstream  sourceFile(inputFileAddress, std::ios::binary);
    std::ofstream  copiedFile(this->GetDirectoryPath() + "_inputParametersCopy.txt",   std::ios::binary);
    copiedFile << sourceFile.rdbuf();
    copiedFile.close();
    sourceFile.close();
    //DictatFiles
    for (PopInt neuronPop : std::ranges::views::iota(0,neurons->GetTotalPopulations())){
        std::shared_ptr<DictatNeuronPop> dictatPtr = std::dynamic_pointer_cast<DictatNeuronPop>(neurons->GetPop(neuronPop));
        if (dictatPtr!=nullptr){
            dictatPtr->CloseInputStream();
            std::ifstream  sourceFile(dictatPtr->GetInputFileAddress(), std::ios::binary);
            std::ofstream  copiedFile(this->directoryPath + nonIterateTitle + "_DictatNeuronPop_"+std::to_string(dictatPtr->GetId())+"_spikers.txt",   std::ios::binary);
            copiedFile << sourceFile.rdbuf();
            copiedFile.close();
            sourceFile.close();
        }
    }
}

//Advanced Recorder functions

void Recorder::SetNoRecordedNeuronPotentials(const std::vector<std::string>& parameterValues){
    PopInt totalNeuronPops = neurons->GetTotalPopulations();
    for(PopInt recordedPop : std::ranges::views::iota(0, std::min(totalNeuronPops, static_cast<PopInt>(parameterValues.size())))){
        neuronPotentialsToRecord.at(recordedPop) = std::stol(parameterValues.at(recordedPop));
        if((neuronPotentialsToRecord.at(recordedPop) >= neurons->GetNeuronsPop(recordedPop)) || (neuronPotentialsToRecord.at(recordedPop) < 0)){
            std::cout << "Potentials: Tracking all neurons of population " << recordedPop << "\n";
            neuronPotentialsToRecord.at(recordedPop) = neurons->GetNeuronsPop(recordedPop);
        }
    }
}

void Recorder::SetNoCurrentContribution(const std::vector<std::string>& parameterValues) {
	PopInt totalNeuronPops = neurons->GetTotalPopulations();
	for (PopInt targetPop : std::ranges::views::iota(0, std::min(totalNeuronPops, static_cast<PopInt>(parameterValues.size())))) {
		currentContributionsToRecord.at(targetPop) = std::min(neurons->GetNeuronsPop(targetPop), std::stol(parameterValues.at(targetPop)));
		if (currentContributionsToRecord.at(targetPop) < 0) {
			std::cout << "Current contribution of all neurons in population " << targetPop << "\n";
			currentContributionsToRecord.at(targetPop) = neurons->GetNeuronsPop(targetPop);
		}
	}
	if (static_cast<PopInt>(parameterValues.size()) > totalNeuronPops) {
		if (parameterValues.at(totalNeuronPops)[0] != '#')// The comments in the code are not meant to be read upon execution
			initialCurrent = static_cast<int>(std::round(std::stod(parameterValues.at(totalNeuronPops)) / infoGlobal->dtTimestep));
	}
}

void Recorder::SetNoRasterplotNeurons(const std::vector<std::string>& parameterValues){
    PopInt totalNeuronPops = neurons->GetTotalPopulations();
    for(PopInt neuronPop : std::ranges::views::iota(0,std::min(totalNeuronPops,static_cast<PopInt>(parameterValues.size())))){
        noRasterPlotNeurons.at(neuronPop) = std::stol(parameterValues.at(neuronPop));
        if((noRasterPlotNeurons.at(neuronPop)  > neurons->GetNeuronsPop(neuronPop)) ||
           (noRasterPlotNeurons.at(neuronPop) < 0)){
            std::cout << "Rasterplot: Tracking all neurons of population "<< neuronPop << "\n";
            noRasterPlotNeurons.at(neuronPop)  = neurons->GetNeuronsPop(neuronPop);
        }
    }
	if (static_cast<PopInt>(parameterValues.size()) > totalNeuronPops) {
		if (parameterValues.at(totalNeuronPops)[0] != '#')
			startRecordingTime = static_cast<int>(std::round(std::stod(parameterValues.at(totalNeuronPops)) / infoGlobal->dtTimestep));
	}
}


void Recorder::SetNoRecordedHeteroSynapticProfilesPerTrackedNeuronPerPop(const std::vector<std::string>& parameterValues) {
    PopInt totalNeuronPops = neurons->GetTotalPopulations();
    for(PopInt neuronPop : std::ranges::views::iota(0,totalNeuronPops)){//min() only makes sense if you remove the hash
        heteroSynTracker.at(neuronPop).first = std::stol(parameterValues.at(2*neuronPop));
        heteroSynTracker.at(neuronPop).second = std::stol(parameterValues.at(2*neuronPop+1));
        // if((heteroSynTracker.at(neuronPop).second  > neurons->GetPop(neuronPop)->GetNoSynapses()) || (heteroSynTracker.at(neuronPop).second < 0)){
        //     heteroSynTracker.at(neuronPop).second=neurons->GetPop(neuronPop)->GetNoSynapses();
        // }
        if((heteroSynTracker.at(neuronPop).first  > neurons->GetNeuronsPop(neuronPop)) || (heteroSynTracker.at(neuronPop).first < 0)){
            heteroSynTracker.at(neuronPop).first=neurons->GetNeuronsPop(neuronPop);
        }
    }
    if (parameterValues.size() > (2 * static_cast<size_t>(totalNeuronPops))) {//This requires an extra number to record
        // std::cout << parameterValues.size() << "_" << 2 * static_cast<size_t>(totalNeuronPops);
        this->heteroRecordPerSteps = std::stol(parameterValues.at(2 * static_cast<size_t>(totalNeuronPops)));
    } //else heteroRecordingPerSteps is 1
}


void Recorder::WriteDataHeaderHeatmap() {
	if (recordedHeatmap == 0)
		return;
	PopInt totalNeuronPops = this->neurons->GetTotalPopulations();
	int dim = infoGlobal->dimensions;
	if (dim != 1 && dim != 2) {
		std::cout << "Heatmap can only be used in 1D or 2D\n";
		throw "Heatmap can only be used in 1D or 2D\n";
	}


	for (PopInt neuronPop : std::ranges::views::iota(0,totalNeuronPops)) {
        std::ofstream heatmapFile;
        AllocateAndAssignStreamBuffer(heatmapFile);
        heatmapFile.open(this->GetHeatmapFilename() + std::to_string(neuronPop)+".dat", std::ofstream::out | std::ofstream::trunc);
        fileStreams.heatmapStreamVector.push_back(std::move(heatmapFile));
		WriteHeader(fileStreams.heatmapStreamVector.back());
		fileStreams.heatmapStreamVector.back() << "# Population:" + std::to_string(neuronPop) + "\n";
		fileStreams.heatmapStreamVector.back() << "# Dimension:" + std::to_string(dim) + "\n";
		fileStreams.heatmapStreamVector.back() << "# L=" + std::to_string(infoGlobal->xAxisLength) + "\n";
		fileStreams.heatmapStreamVector.back() << "#************************************\n";
		if (dim == 1) {
			fileStreams.heatmapStreamVector.back() << "# column 1 : t (secs.) \n";
			fileStreams.heatmapStreamVector.back() << "# columnn n+1: r_" + std::to_string(neuronPop) + "(Hz) for neurons in position:\n";
			fileStreams.heatmapStreamVector.back() << "#\t\t (n-1) * " + std::to_string(infoGlobal->xAxisLength / recordedHeatmap) + "  <  x  <  n * " + std::to_string(infoGlobal->xAxisLength / recordedHeatmap) + "\n";
		}
		else {
			fileStreams.heatmapStreamVector.back() << "# column 1 : t (secs.) \n";
			fileStreams.heatmapStreamVector.back() << "# columnn n+1: r_" + std::to_string(neuronPop) + "(Hz) for neurons in position:\n";
			fileStreams.heatmapStreamVector.back() << "#\t\t (n-1)%" + std::to_string(recordedHeatmap) + " * " + std::to_string(infoGlobal->xAxisLength / recordedHeatmap) + "  <  x  <  n%" + std::to_string(recordedHeatmap) + " * " + std::to_string(infoGlobal->xAxisLength / recordedHeatmap) + "\n";
			fileStreams.heatmapStreamVector.back() << "#\t\t (n-1)//" + std::to_string(recordedHeatmap) + " * " + std::to_string(infoGlobal->xAxisLength / recordedHeatmap) + "  <  y  <  n//" + std::to_string(recordedHeatmap) + " * " + std::to_string(infoGlobal->xAxisLength / recordedHeatmap) + "\n";
		}
		//legend header
		if (dim == 1) {
			fileStreams.heatmapStreamVector.back() << "t\t";
			for (int xPosition : std::ranges::views::iota(0,recordedHeatmap)) {
				fileStreams.heatmapStreamVector.back() <<"["<< static_cast<double>(xPosition) / static_cast<double>(recordedHeatmap)*infoGlobal->xAxisLength<<"]\t";
			}
		}
		else if (dim == 2) {
			fileStreams.heatmapStreamVector.back() << "t\t";
			for (int xPosition : std::ranges::views::iota(0,recordedHeatmap)) {
				for (int yPosition : std::ranges::views::iota(0,recordedHeatmap)) {
					fileStreams.heatmapStreamVector.back() << "[" << static_cast<double>(xPosition) / static_cast<double>(recordedHeatmap)*infoGlobal->xAxisLength <<","<< static_cast<double>(yPosition) / static_cast<double>(recordedHeatmap)*infoGlobal->xAxisLength << "]\t";
				}
			}
		}
		fileStreams.heatmapStreamVector.back() << "\n#************************************\n";
		for (int kIndex : std::ranges::views::iota (1, static_cast<int>(2 + pow(recordedHeatmap, dim)))){
			fileStreams.heatmapStreamVector.back() << "# " + std::to_string(kIndex) << "  \t";
        }
		fileStreams.heatmapStreamVector.back() << std::endl;
	}
}

void Recorder::WriteDataHeaderAverages(){

    PopInt totalNeuronPops = this->neurons->GetTotalPopulations();
    int columnTotal = 1;

    //I have not reformatted this function yet because .
    AllocateAndAssignStreamBuffer(this->fileStreams.averagesFileStream);
    this->fileStreams.averagesFileStream.open (this->GetDataFilename(), std::ofstream::out | std::ofstream::trunc);
    WriteHeader(this->fileStreams.averagesFileStream);
    this->fileStreams.averagesFileStream << "#************************************\n";
    this->fileStreams.averagesFileStream << "#Columns: \n";
    this->fileStreams.averagesFileStream << "#" +std::to_string(columnTotal)+ " t (secs.) \n"; 
    columnTotal++;
    for(PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
        this->fileStreams.averagesFileStream << "#" +std::to_string(columnTotal)+ " V_" << neuronPop << " (mV)\n";
        columnTotal++;
    }

    for(PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
        this->fileStreams.averagesFileStream << "#" +std::to_string(columnTotal)+ " r_" << neuronPop << " (Hz)\n"; 
        columnTotal++;
    }

    for(PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
        this->fileStreams.averagesFileStream << "#" +std::to_string(columnTotal)+ " mu_Ext_" << neuronPop << "/tau_m (dmV/sec) \n"; 
        columnTotal++;
    }

    for(PopInt targetPop : std::ranges::views::iota(0, totalNeuronPops)){
        for(PopInt sourcePop : std::ranges::views::iota(0, totalNeuronPops)){
            this->fileStreams.averagesFileStream << "# " +std::to_string(columnTotal)+ " mu_" << std::to_string(targetPop) << "_" << std::to_string(sourcePop) << "/tau_m (dmV/sec) \n"; 
            columnTotal++;
        }
    }
    for(PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
        this->fileStreams.averagesFileStream << "#" +std::to_string(columnTotal)+ " mu_" << neuronPop << "/tau_m (dmV/sec) \n"; 
        columnTotal++;
    }
    for(PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
        this->fileStreams.averagesFileStream << "#" +std::to_string(columnTotal)+ " Quenched fluctuations of mu_" << neuronPop << "/tau_m (dmV/sec) = sqrt(PopAver(TempAver(mu_i)^2) - PopAver(TempAver(mu_i))^2) \n"; 
        columnTotal++;
    }

    for(PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
        this->fileStreams.averagesFileStream << "#" +std::to_string(columnTotal)+ " Temporal fluctuations of mu_" << neuronPop << "/tau_m*sqrt(dt) (dmV/sqrt(sec)) = sqrt(dt)*sqrt(PopAver(TempAver(mu_i^2) - TempAver(mu_i)^2))\n"; 
        columnTotal++;
    }

	this->fileStreams.averagesFileStream << "t\t";
	for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
		this->fileStreams.averagesFileStream << "V_" << neuronPop<<"\t";
    }
	for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
		this->fileStreams.averagesFileStream << "r_" << neuronPop << "\t";
    }
	for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
		this->fileStreams.averagesFileStream << "mu_ext_" << neuronPop << "\t";
    }
	for (PopInt targetPop : std::ranges::views::iota(0, totalNeuronPops)) {
		for (PopInt sourcePop : std::ranges::views::iota(0, totalNeuronPops)){
			this->fileStreams.averagesFileStream << "mu_" << targetPop << "_"<< sourcePop << "\t";
        }
	}
	for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
		this->fileStreams.averagesFileStream << "mu_tot_" << neuronPop << "\t";
    }
	for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
        this->fileStreams.averagesFileStream << "sigma_quenched_" << neuronPop << "\t";
    }
	for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
		this->fileStreams.averagesFileStream << "sigma_t_" << neuronPop << "\t";
    }
	this->fileStreams.averagesFileStream << "\n#************************************\n";

    for(int column : std::ranges::views::iota (1, columnTotal)){
        this->fileStreams.averagesFileStream << "#"+std::to_string(column) << "\t\t";
    }
    this->fileStreams.averagesFileStream << std::endl;
}

void Recorder::WriteDataHeaderRasterplot(){

    if(!recordRasterPlot){
        return;
    }
    AllocateAndAssignStreamBuffer(this->fileStreams.rasterplotFileStream);
    this->fileStreams.rasterplotFileStream.open (GetRasterplotFilename(), std::ofstream::out | std::ofstream::trunc);
    WriteHeader(this->fileStreams.rasterplotFileStream);
    this->fileStreams.rasterplotFileStream << "#1 t (secs.) \t #2 neuron_id \t #3 neuron_pop_id \t\n";
	this->fileStreams.rasterplotFileStream << "# Note that the recorded neurons are equidistant within the population. Therefore, neuron_id are not necessarily successive numbers\n";
	this->fileStreams.rasterplotFileStream << "Spike_t\tNeuron_id\tPop_id\n";
	this->fileStreams.rasterplotFileStream << "#************************************"<<std::endl;
}

void Recorder::WriteDataHeaderPotential(){

    if(!recordNeuronPotentials){
        return;
    }

    PopInt  totalNeuronPops = neurons->GetTotalPopulations();
    AllocateAndAssignStreamBuffer(this->fileStreams.potentialFileStream);
    this->fileStreams.potentialFileStream.open(GetPotentialFilename(), std::ofstream::out | std::ofstream::trunc);

    WriteHeader(this->fileStreams.potentialFileStream);
    this->fileStreams.potentialFileStream << "#1 t (secs.)\t 2-"<<1+std::reduce(neuronPotentialsToRecord.begin(), neuronPotentialsToRecord.end())<<" V_pop_id (mV) \n";

	this->fileStreams.potentialFileStream << "t\t";
    for(PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
        for(PopInt neuron : std::ranges::views::iota(0,neuronPotentialsToRecord.at(neuronPop))){
            this->fileStreams.potentialFileStream << "V_" << neuronPop << "_" << neuron <<  "\t";
        }
    }
	this->fileStreams.potentialFileStream << "\n#************************************"<<std::endl;
}

void Recorder::WriteDataHeaderCurrents(){
    if(!recordNeuronPotentials){
        return;
    }


    PopInt totalNeuronPops = neurons->GetTotalPopulations();
    AllocateAndAssignStreamBuffer(this->fileStreams.currentsFileStream);
    this->fileStreams.currentsFileStream.open(GetCurrentsFilename(), std::ofstream::out | std::ofstream::trunc);

    WriteHeader(this->fileStreams.currentsFileStream);
    this->fileStreams.currentsFileStream << "#1 t (sec)\t 2-"<<1+ std::reduce(neuronPotentialsToRecord.begin(), neuronPotentialsToRecord.end())<<" mu_pop_id / tau_m(dmV / sec)\n";

	this->fileStreams.currentsFileStream << "t\t";
	for(PopInt neuronPop : std::ranges::views::iota(0,totalNeuronPops)){
        for(NeuronInt neuron : std::ranges::views::iota(0,neuronPotentialsToRecord.at(neuronPop))){
            this->fileStreams.currentsFileStream << "mu_" << neuronPop << "_" << neuron << "\t";
            // std::cout << neuron<<"_"<<neuronPop << std::endl;
        }
    }
	this->fileStreams.currentsFileStream << "\n#************************************"<<std::endl;

}

void Recorder::WriteDataHeaderCurrentsContribution() {
    if(!recordCurrentContributions){
		return;
    }
	PopInt totalNeuronPops = neurons->GetTotalPopulations();
	long HeaderIndex{};
	NeuronInt neuronNumber{};
    AllocateAndAssignStreamBuffer(this->fileStreams.cCurrentsFileStream);
	this->fileStreams.cCurrentsFileStream.open(GetCurrentCrontributionFilename(), std::ofstream::out | std::ofstream::trunc);

	WriteHeader(this->fileStreams.cCurrentsFileStream);
	this->fileStreams.cCurrentsFileStream << "# 1 t (sec)\n";
	HeaderIndex = 2;
	for (PopInt targetPop : std::ranges::views::iota (0, totalNeuronPops)) {
		neuronNumber = currentContributionsToRecord.at(targetPop);
		if (neuronNumber > 0) {
			for (PopInt sourcePop : std::ranges::views::iota(0, totalNeuronPops)) {
				this->fileStreams.cCurrentsFileStream <<"# "<< HeaderIndex <<"-"<< HeaderIndex +neuronNumber -1 <<" Irecurrent_pop"<< targetPop<<"<-neuronPop"<< sourcePop<<"_neuronid / tau_m (mV/s)\n";
				HeaderIndex = HeaderIndex + neuronNumber;
			}
			this->fileStreams.cCurrentsFileStream << "# "<<HeaderIndex << "-" << HeaderIndex + neuronNumber -1 << " Ifeedforward_pop"<< targetPop <<"_id / tau_m (mV/s)\n";
			HeaderIndex = HeaderIndex + neuronNumber;
		}
	}
	this->fileStreams.cCurrentsFileStream << "# Note that the recorded neurons are equidistant within the population. Therefore, neuron_id are not necessarily successive numbers\n";

	this->fileStreams.cCurrentsFileStream << "t\t";
	for (PopInt targetPop : std::ranges::views::iota (0, totalNeuronPops)) {
		for (PopInt sourcePop : std::ranges::views::iota(0, totalNeuronPops)) {
			for (NeuronInt recordedNeuron : std::ranges::views::iota(0,currentContributionsToRecord.at(targetPop))) {
                //This is an integer type division, it does not need floor()
				this->fileStreams.cCurrentsFileStream << "IRec" << targetPop << "_" << sourcePop << "_" << recordedNeuron*neurons->GetNeuronsPop(targetPop) / currentContributionsToRecord.at(targetPop) << "\t";
			}
		}
		for (NeuronInt recordedNeuron : std::ranges::views::iota(0,currentContributionsToRecord.at(targetPop))) {
            //This is an integer type division, it does not need floor() 
			this->fileStreams.cCurrentsFileStream << "Iffd" << targetPop << "_" << recordedNeuron*neurons->GetNeuronsPop(targetPop) / currentContributionsToRecord.at(targetPop) << "\t";
		}
	}
	this->fileStreams.cCurrentsFileStream << "\n#************************************"<<std::endl;
}


void Recorder::WriteDataHeaderSynapseStates() {
    if(!trackSynapses){
        return;
    }
    AllocateAndAssignStreamBuffer(this->fileStreams.synStatesFileStream);
    this->fileStreams.synStatesFileStream.open (this->GetSynapseStateFilename(), std::ofstream::out | std::ofstream::trunc);
    WriteHeader(this->fileStreams.synStatesFileStream);
    this->fileStreams.synStatesFileStream << "#************************************\n";
    this->fileStreams.synStatesFileStream << "#Columns: \n";
    this->fileStreams.synStatesFileStream << "#1 t (secs.) \n";
	this->fileStreams.synStatesFileStream << synapses->GetDataHeader(2);
	this->fileStreams.synStatesFileStream << "t\t"<<synapses->GetUnhashedDataHeader()<< "\n";
    this->fileStreams.synStatesFileStream << "# Attention, Synaptic state data is spike-induced : at each time step, only synapses from which the presynaptic neuron has spiked are measured. Data only tested for CurrentSynapse and MongilloSynapse\n";
    this->fileStreams.synStatesFileStream << "#************************************"<<std::endl;
}


void Recorder::WriteDataHeaderHeteroSynapses(){

    if(!recordHeteroSynapses){
        return;
    }

    PopInt totalNeuronPops = neurons->GetTotalPopulations();
    AllocateAndAssignStreamBuffer(this->fileStreams.heteroSynFileStream);
    this->fileStreams.heteroSynFileStream.open(GetHeteroSynapseStateFilename(), std::ofstream::out | std::ofstream::trunc);
    WriteHeader(this->fileStreams.heteroSynFileStream);

    for(PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
        BindNoHeteroSynapsesPerPop(neuronPop);
        if (heteroSynTracker.at(neuronPop).first != 0 && heteroSynTracker.at(neuronPop).second != 0 && this->neurons->GetPop(neuronPop)->HasPlasticityModel()){
            this->fileStreams.heteroSynFileStream << "#Pop. "<< neuronPop << " profile -> "<<this->neurons->GetPop(neuronPop)->GetIndividualSynapticProfileHeaderInfo() <<" \n";
        }
    }

    this->fileStreams.heteroSynFileStream << "\n#************************************\n";
    this->fileStreams.heteroSynFileStream << "#1 t (secs.)\t 2-"<<1+std::accumulate(heteroSynTracker.begin(), heteroSynTracker.end(), 0, [](int accumulator,std::pair<NeuronInt, signed long> pop){return accumulator + (pop.first*pop.second);})<<" Profile_pop_id_neuron_id \n";
    this->fileStreams.heteroSynFileStream << "t\t";

    for(PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
        if ( heteroSynTracker.at(neuronPop).second == 0 ||heteroSynTracker.at(neuronPop).first ==0|| !this->neurons->GetPop(neuronPop)->HasPlasticityModel()) {
            continue;
        }
        for(NeuronInt neuron : std::ranges::views::iota(0,heteroSynTracker.at(neuronPop).first)) {
            for (signed long synapse : std::ranges::views::iota(0,heteroSynTracker.at(neuronPop).second)) {
                this->fileStreams.heteroSynFileStream << "Profile_" << neuronPop << "_" << neuron << "_" << synapse <<  "\t";
            }
        }
    }
    this->fileStreams.heteroSynFileStream << "\n#************************************"<<std::endl;
}

void Recorder::WriteDataHeaderHeteroSynapsesOverall(){

    if(!recordHeteroSynapses){
        return;
    }

    //std::cout << "The file has been properly created!!!!\n";
    PopInt totalNeuronPops = neurons->GetTotalPopulations();
    AllocateAndAssignStreamBuffer(this->fileStreams.hSOverallFileStream);
    this->fileStreams.hSOverallFileStream.open(GetOverallHeteroSynapseStateFilename(), std::ofstream::out | std::ofstream::trunc);

    WriteHeader(this->fileStreams.hSOverallFileStream);
    for(PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
        BindNoHeteroSynapsesPerPop(neuronPop);
        if (heteroSynTracker.at(neuronPop).first != 0 &&heteroSynTracker.at(neuronPop).second != 0 && this->neurons->GetPop(neuronPop)->HasPlasticityModel()){
            this->fileStreams.hSOverallFileStream << "#Pop. "<< neuronPop << " Overall Profile -> "<<this->neurons->GetPop(neuronPop)->GetOverallSynapticProfileHeaderInfo() <<" \n";
        }
    }
    this->fileStreams.hSOverallFileStream << "\n#************************************\n";

    this->fileStreams.hSOverallFileStream << "#1 t (secs.)\t 2-"<<1+std::accumulate(heteroSynTracker.begin(), heteroSynTracker.end(), 0, [](int accumulator,std::pair<NeuronInt, signed long> pop){return accumulator + (pop.first*pop.second);})<<" Profile_pop_id_neuron_id \n";

    this->fileStreams.hSOverallFileStream << "t\t";

    for(PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
        if ( heteroSynTracker.at(neuronPop).second == 0 ||heteroSynTracker.at(neuronPop).first ==0|| !this->neurons->GetPop(neuronPop)->HasPlasticityModel()) {
            continue;
        }
        for(NeuronInt neuron : std::ranges::views::iota(0,heteroSynTracker.at(neuronPop).first)) {
            this->fileStreams.hSOverallFileStream << "OverallProfile_" << neuronPop << "_" << neuron <<  "\t";
        }
    }

    this->fileStreams.hSOverallFileStream << "\n#************************************"<<std::endl;
}


void Recorder::WriteDataHeader(){
    WriteDataHeaderRasterplot();
    WriteDataHeaderAverages();
    WriteDataHeaderSynapseStates();
    WriteDataHeaderPotential();
    WriteDataHeaderCurrents();
	WriteDataHeaderHeatmap();
	WriteDataHeaderCurrentsContribution();
	WriteDataHeaderHeteroSynapses();
	WriteDataHeaderHeteroSynapsesOverall();
    ResetStatistics();
}

void Recorder::ResetStatistics()
{
    PopInt totalNeuronPops = this->neurons->GetTotalPopulations();

    std::fill(currentBin.potential.begin(), currentBin.potential.end(), 0.0);
    std::fill(currentBin.spikerRatio.begin(), currentBin.spikerRatio.end(), 0.0);
    std::fill(currentBin.externalCurrent.begin(), currentBin.externalCurrent.end(), 0.0);
    std::fill(currentBin.totalCurrentSquared_mean.begin(), currentBin.totalCurrentSquared_mean.end(), 0.0);

    for(PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
        std::fill(currentBin.neuronTotalCurrentMean.at(neuronPop).begin(), currentBin.neuronTotalCurrentMean.at(neuronPop).end(), 0.0);
        std::fill(currentBin.noRecordedSynapses.at(neuronPop).begin(),currentBin.noRecordedSynapses.at(neuronPop).end(), 0);
        std::fill(currentBin.synapticCurrents.at(neuronPop).begin(), currentBin.synapticCurrents.at(neuronPop).end(), 0.0);
        for(PopInt sourcePop : std::ranges::views::iota(0, totalNeuronPops)) {
            //currentBin.noRecordedSynapses.at(neuronPop)[neuronPop2] = 0;
            //currentBin.synapticState.at(neuronPop).at(sourcePop)       = 0.0;
            std::fill(currentBin.synapticState.at(neuronPop).at(sourcePop).begin(), currentBin.synapticState.at(neuronPop).at(sourcePop).end(), 0.0);
        }
		if (recordedHeatmap != 0) {
			std::fill(currentBin.heatmap.at(neuronPop).begin(),currentBin.heatmap.at(neuronPop).end(), 0.0);
		}
    }
}

void Recorder::RecordRasterplot(){

	if (!recordRasterPlot || infoGlobal->timeStep < startRecordingTime){
        return;
    }

    double           dtTimestep = infoGlobal->dtTimestep;
    // int				first_id;
	NeuronInt		neuronsToRecord;//number of neurons to record
	NeuronInt		totalNeurons;//number of neurons in population

    for(PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())){
		neuronsToRecord = noRasterPlotNeurons.at(neuronPop);
		if (neuronsToRecord == 0){
			continue;
        }
        const std::vector<NeuronInt>& spikerVector = neurons->GetSpikers(neuronPop);
		totalNeurons = neurons->GetNeuronsPop(neuronPop);

        for(const NeuronInt spikerID : spikerVector){
            if (spikerID == std::lround((static_cast<double>(totalNeurons) / neuronsToRecord) * std::lround(static_cast<double>(neuronsToRecord * spikerID) / totalNeurons))){//WTH is this
                this->fileStreams.rasterplotFileStream << static_cast<double>(infoGlobal->timeStep)*dtTimestep <<"\t" << spikerID << "\t" << neuronPop << "\n";
            }
        }
    }
}

void Recorder::RecordCurrents(const std::vector<std::vector<double>>& synaptic_dV){

    if(!recordNeuronPotentials){
        return;
    }

    double           dtTimestep = infoGlobal->dtTimestep;
    double           time = static_cast<double>(infoGlobal->timeStep)*dtTimestep;
    PopInt           totalNeuronPops = neurons->GetTotalPopulations();

    SaveDoubleFile(this->fileStreams.currentsFileStream,time,5);
    for(PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)){
        for(NeuronInt recordedNeuron : std::ranges::views::iota(0, neuronPotentialsToRecord.at(neuronPop))){
            SaveDoubleFile(this->fileStreams.currentsFileStream,synaptic_dV.at(neuronPop).at(recordedNeuron)/dtTimestep,5);
        }
    }
    this->fileStreams.currentsFileStream << "\n";
}

void Recorder::RecordPotential(){

    if(!recordNeuronPotentials){
        return;
    }

    double           dtTimestep = infoGlobal->dtTimestep;
    double           time = static_cast<double>(infoGlobal->timeStep)*dtTimestep;
    PopInt           totalNeuronPops = neurons->GetTotalPopulations();

    SaveDoubleFile(this->fileStreams.potentialFileStream,time,5);
    for(PopInt neuronPop : std::ranges::views::iota(0,totalNeuronPops)){
        for(NeuronInt recordedNeuron : std::ranges::views::iota(0, neuronPotentialsToRecord.at(neuronPop))){
            SaveDoubleFile(this->fileStreams.potentialFileStream,neurons->GetPotential(neuronPop,recordedNeuron),5);
        }
    }
    this->fileStreams.potentialFileStream << "\n";
}

void Recorder::RecordCurrentContributions(const std::vector<std::vector<double>>& synaptic_dV) {

	if (!recordCurrentContributions || infoGlobal->timeStep < initialCurrent){
		return;
    }

	PopInt totalNeuronPops = neurons->GetTotalPopulations();
	NeuronInt indexMemory{};


	for (PopInt targetPop : std::ranges::views::iota(0, totalNeuronPops)) {
		for (PopInt sourcePop : std::ranges::views::iota(0, totalNeuronPops)) {
			for (NeuronInt recordedNeuron : std::ranges::views::iota (0,currentContributionsToRecord.at(targetPop))) {
                if(synapses->GetConnectedState(targetPop, sourcePop)){
                    currentContrBin.at(static_cast<size_t>(recordedNeuron) + static_cast<size_t>(indexMemory)) += synapses->GetRecurrentInput(targetPop, sourcePop, (recordedNeuron*neurons->GetNeuronsPop(targetPop) / currentContributionsToRecord.at(targetPop)));
                }
			}
			indexMemory += currentContributionsToRecord.at(targetPop);
		}
		for (NeuronInt recordedNeuron : std::ranges::views::iota (0,currentContributionsToRecord.at(targetPop))) {
			currentContrBin.at(static_cast<size_t>(recordedNeuron) + static_cast<size_t>(indexMemory)) += stimulus->GetSignalMatrixPoint(targetPop, (recordedNeuron*neurons->GetNeuronsPop(targetPop) / currentContributionsToRecord.at(targetPop)));
		}
		indexMemory += currentContributionsToRecord.at(targetPop);
	}

	if ((infoGlobal->timeStep) % this->timeStepsPerBin == 0) {
		SaveDoubleFile(this->fileStreams.cCurrentsFileStream, static_cast<double>(infoGlobal->timeStep)*infoGlobal->dtTimestep, 5);
		for (double storedValue : currentContrBin) {
			SaveDoubleFile(this->fileStreams.cCurrentsFileStream, storedValue, 5);
		}
        std::fill(currentContrBin.begin(), currentContrBin.end(), 0.0);
		this->fileStreams.cCurrentsFileStream << "\n";
	}
}

void Recorder::RecordSynapseStates(){

    PopInt totalNeuronPops{neurons->GetTotalPopulations()};

    if(!trackSynapses){
        return;
    }

    if((infoGlobal->timeStep+1)%this->timeStepsPerBin == 0) {
        SaveDoubleFile(this->fileStreams.synStatesFileStream,static_cast<double>(infoGlobal->timeStep)*infoGlobal->dtTimestep,4);
        //Record synaptic strengths
        for(PopInt sourcePop : std::ranges::views::iota(0, totalNeuronPops)) {
            for(PopInt targetPop : std::ranges::views::iota(0, totalNeuronPops)) {
                if(currentBin.noRecordedSynapses.at(targetPop).at(sourcePop) > 0) {
                    for (double recordedValue : currentBin.synapticState.at(targetPop).at(sourcePop)) {
                        SaveDoubleFile(this->fileStreams.synStatesFileStream,recordedValue /static_cast<double>(currentBin.noRecordedSynapses.at(targetPop).at(sourcePop)),6);
                    }
                } else {
                    for(signed dataLength : std::ranges::views::iota(0, static_cast<signed>(currentBin.synapticState.at(targetPop).at(sourcePop).size()))){
                        (void)dataLength;
                        SaveDoubleFile(this->fileStreams.synStatesFileStream,NAN,6);
                        // SaveDoubleFile(this->fileStreams.synStatesFileStream,0,6);
                    }
                }
            }
        }
        this->fileStreams.synStatesFileStream << "\n";
    }
}

void Recorder::RecordAverages(){

    PopInt    noTotalPopulations       = neurons->GetTotalPopulations();
    double timeBinSize {static_cast<double>(this->timeStepsPerBin)*infoGlobal->dtTimestep};
    double binSizeInSteps      {static_cast<double>(timeStepsPerBin)};
    double noNeurons;

    std::vector<double> means(noTotalPopulations, 0.0);
    std::vector<double> popAveraged_SquaredMeanTime(noTotalPopulations, 0.0);

    for(PopInt neuronPop : std::ranges::views::iota(0,noTotalPopulations)){
        noNeurons = neurons->GetNeuronsPop(neuronPop);
        for(NeuronInt neuron : std::ranges::views::iota (0, neurons->GetNeuronsPop(neuronPop))){
            means.at(neuronPop)                  += currentBin.neuronTotalCurrentMean.at(neuronPop).at(neuron)/(noNeurons*binSizeInSteps);
            popAveraged_SquaredMeanTime.at(neuronPop)+= pow(currentBin.neuronTotalCurrentMean.at(neuronPop).at(neuron)/binSizeInSteps,2)/noNeurons; // population average over square of temporal average [(I_alpha^i)^2]
        }
    }


    if((infoGlobal->timeStep)%this->timeStepsPerBin == 0) {
        SaveDoubleFile(this->fileStreams.averagesFileStream,static_cast<double>(infoGlobal->timeStep)*infoGlobal->dtTimestep,4);
        //Record average potential
        for(double popRecord : currentBin.potential) {
            SaveDoubleFile(this->fileStreams.averagesFileStream,popRecord/binSizeInSteps,3);
        }
        //Record firing rate
        for(double popRecord : currentBin.spikerRatio) {
            SaveDoubleFile(this->fileStreams.averagesFileStream,popRecord/timeBinSize,3);
        }
        //Record external input 'current' (in mV/sec)
        for(double popRecord : currentBin.externalCurrent) {
            SaveDoubleFile(this->fileStreams.averagesFileStream,popRecord/binSizeInSteps,3);
        }
        //Record population averaged synaptic input currents mu_ij (in mV/sec)
        for(std::vector<double>& targetPopRecord : currentBin.synapticCurrents) {
            for(double record : targetPopRecord) {
                SaveDoubleFile(this->fileStreams.averagesFileStream, record/binSizeInSteps,3);
            }
        }
        //Record population averaged total synaptic input currents mu_i (in mV/sec)
        //for(int i = 0; i < noTotalPopulations; i++)
        //{
        //    value   = currentBin.totalCurrent_mean[i]/n_aver;
        //    SaveDoubleFile(this->FileStreams.averagesFileStream,value,3);
        //}

        //Record check of mean total input currents
        for(double popRecord : means) {
            SaveDoubleFile(this->fileStreams.averagesFileStream,popRecord,3);
        }

        //Quenched fluctuations (due to varying input over population)
        for(PopInt neuronPop : std::ranges::views::iota (0, noTotalPopulations)) {
            SaveDoubleFile(this->fileStreams.averagesFileStream,sqrt(popAveraged_SquaredMeanTime.at(neuronPop)-pow(means.at(neuronPop),2)),3);
        }

        //Temporal fluctuations
        for(PopInt neuronPop : std::ranges::views::iota (0, noTotalPopulations)) {
            SaveDoubleFile(this->fileStreams.averagesFileStream,sqrt(currentBin.totalCurrentSquared_mean.at(neuronPop)/binSizeInSteps - popAveraged_SquaredMeanTime.at(neuronPop))*infoGlobal->dtSqrt,3);
        }

        this->fileStreams.averagesFileStream << "\n";
        ResetStatistics();
    }
}
//
void Recorder::RecordHeatmap() {

	if (recordedHeatmap == 0){
		return;
    }

	PopInt    noTotalPopulations = neurons->GetTotalPopulations();
	double dtTimestep = infoGlobal->dtTimestep;
	double timeBinSize = static_cast<double>(this->timeStepsPerBin)*dtTimestep;
	// double n_aver = static_cast<double>(timeStepsPerBin);
	// double n;

	if ((infoGlobal->timeStep) % this->timeStepsPerBin == 0){
		for (PopInt neuronPop : std::ranges::views::iota(0,noTotalPopulations)) {
			SaveDoubleFile(fileStreams.heatmapStreamVector.at(neuronPop), static_cast<double>(infoGlobal->timeStep)*dtTimestep, 4);
			for (double pixel : currentBin.heatmap.at(neuronPop)) {
				SaveDoubleFile(fileStreams.heatmapStreamVector.at(neuronPop), pixel / timeBinSize, 3);
			}
			fileStreams.heatmapStreamVector.at(neuronPop) << "\n";
		}
	}
}

void Recorder::RecordHeteroSynapses() {

    if(!recordHeteroSynapses) {
        return;
    }

    SaveDoubleFile(this->fileStreams.heteroSynFileStream,static_cast<double>(infoGlobal->timeStep)*infoGlobal->dtTimestep,5);

    for(PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())){
        if (!this->neurons->GetPop(neuronPop)->HasPlasticityModel() || heteroSynTracker.at(neuronPop).second == 0|| heteroSynTracker.at(neuronPop).first == 0) {
            continue;
        }
        for(NeuronInt neuron : std::ranges::views::iota(0,heteroSynTracker.at(neuronPop).first)) {
            for (signed long synapse : std::ranges::views::iota(0,heteroSynTracker.at(neuronPop).second)) {
                SaveTupleOfDoublesFile(this->fileStreams.heteroSynFileStream, this->neurons->GetPop(neuronPop)->GetIndividualSynapticProfile(neuron, synapse), 5);
            }
        }
    }
    this->fileStreams.heteroSynFileStream << "\n";
}

void Recorder::RecordHeteroSynapsesOverall() {
    if(!recordHeteroSynapses) {
        return;
    }

    SaveDoubleFile(this->fileStreams.hSOverallFileStream,static_cast<double>(infoGlobal->timeStep)*infoGlobal->dtTimestep,5);

    for(PopInt neuronPop : std::ranges::views::iota(0, neurons->GetTotalPopulations())){
        if (!this->neurons->GetPop(neuronPop)->HasPlasticityModel() || heteroSynTracker.at(neuronPop).second == 0|| heteroSynTracker.at(neuronPop).first == 0) {
            continue;
        }
        for(NeuronInt neuron : std::ranges::views::iota(0,heteroSynTracker.at(neuronPop).first)) {
            SaveTupleOfDoublesFile(this->fileStreams.hSOverallFileStream, this->neurons->GetPop(neuronPop)->GetOverallSynapticProfile(neuron), 5);
            //Here is selecting only the average weight per neuron, with precision 5 digits.
        }
    }
    this->fileStreams.hSOverallFileStream << "\n";
}

void Recorder::Record(const std::vector<std::vector<double>>& synaptic_dV) {
	PopInt totalNeuronPops = this->neurons->GetTotalPopulations();
	double noNeurons;
	// int Xbin, Ybin;

	for (PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops)) {
		for (NeuronInt neuron : std::ranges::views::iota(0, neurons->GetNeuronsPop(neuronPop))) {
			noNeurons = static_cast<double>(this->neurons->GetNeuronsPop(neuronPop));
			currentBin.potential.at(neuronPop) += this->neurons->GetPotential(neuronPop, neuron) / noNeurons;
			currentBin.externalCurrent.at(neuronPop) += this->stimulus->GetSignalMatrixPoint(neuronPop, neuron) / (noNeurons * infoGlobal->dtTimestep);
			currentBin.totalCurrentSquared_mean.at(neuronPop) += pow(synaptic_dV.at(neuronPop).at(neuron) / infoGlobal->dtTimestep, 2.0) / noNeurons;
			//record per neuron
			currentBin.neuronTotalCurrentMean.at(neuronPop).at(neuron) += synaptic_dV.at(neuronPop).at(neuron) / infoGlobal->dtTimestep;
		}
	}
	for (PopInt sourcePop : std::ranges::views::iota(0,totalNeuronPops)) {
		const std::vector<NeuronInt>& spikers = this->neurons->GetSpikersPrevDt(sourcePop);
		currentBin.spikerRatio.at(sourcePop) += static_cast<double>(spikers.size()) / static_cast<double>(this->neurons->GetNeuronsPop(sourcePop));
		if (recordedHeatmap != 0) {
			if (infoGlobal->dimensions == 2) {
                    std::for_each(PAR_UNSEQ,spikers.begin(), spikers.end(),[this, sourcePop](double spiker){
					double Xbin = static_cast<int>(neurons->GetXPosition(sourcePop, spiker) * recordedHeatmap / infoGlobal->xAxisLength);
					double Ybin = static_cast<int>(neurons->GetYPosition(sourcePop, spiker) * recordedHeatmap / infoGlobal->yAxisLength);
					currentBin.heatmap.at(sourcePop).at(static_cast<size_t>(recordedHeatmap) * Ybin + Xbin) += 1 / static_cast<double>(densimap.at(sourcePop).at(static_cast<size_t>(recordedHeatmap) * Ybin + Xbin));
                });
				// for (NeuronInt spiker : spikers) {
				// 	Xbin = static_cast<int>(neurons->GetXPosition(sourcePop, spiker) * recordedHeatmap / infoGlobal->xAxisLength);
				// 	Ybin = static_cast<int>(neurons->GetYPosition(sourcePop, spiker) * recordedHeatmap / infoGlobal->yAxisLength);
                //     currentBin.heatmap.at(sourcePop).at(static_cast<size_t>(recordedHeatmap) * Ybin + Xbin) += 1 / static_cast<double>(densimap.at(sourcePop).at(static_cast<size_t>(recordedHeatmap) * Ybin + Xbin));
				// }
			} else {
                std::for_each(PAR_UNSEQ,spikers.begin(), spikers.end(),[this, sourcePop](double spiker){
					double Xbin = static_cast<int>(neurons->GetXPosition(sourcePop, spiker) * recordedHeatmap / infoGlobal->xAxisLength);
					currentBin.heatmap.at(sourcePop).at(Xbin) += 1 / static_cast<double> (densimap.at(sourcePop).at(Xbin));
                });
				// for (NeuronInt spiker : spikers) {
				// 	Xbin = static_cast<int>(neurons->GetXPosition(sourcePop, spiker) * recordedHeatmap / infoGlobal->xAxisLength);
				// 	currentBin.heatmap.at(sourcePop).at(Xbin) += 1 / static_cast<double> (densimap.at(sourcePop).at(Xbin));
				// }
			}
		}
		for (PopInt targetPop : std::ranges::views::iota(0, totalNeuronPops)){
            if(synapses->GetConnectedState(targetPop, sourcePop)){
			    currentBin.synapticCurrents.at(targetPop).at(sourcePop) += (this->synapses->GetCumulatedDV(targetPop, sourcePop) / static_cast<double>(this->neurons->GetNeuronsPop(targetPop))) / infoGlobal->dtTimestep;
            }
        }
		if (!trackSynapses){
			continue;
        }
		for (NeuronInt spiker : spikers) {
			for (PopInt targetPop : std::ranges::views::iota(0, totalNeuronPops)) {
                if (synapses->GetConnectedState(targetPop, sourcePop)){
                    std::vector<double> synapticState {this->synapses->GetSynapticState(targetPop, sourcePop, spiker)};
                    std::transform(PAR_UNSEQ,currentBin.synapticState.at(targetPop).at(sourcePop).begin(), currentBin.synapticState.at(targetPop).at(sourcePop).end(),synapticState.begin(), currentBin.synapticState.at(targetPop).at(sourcePop).begin(), std::plus<double>());
                    //currentBin.synapticState.at(targetPop).at(sourcePop) += this->synapses->GetSynapticState(targetPop, sourcePop, spiker);//not valid with vectors 
                    currentBin.noRecordedSynapses.at(targetPop).at(sourcePop) += this->synapses->GetNoTargetedNeurons(targetPop, sourcePop, spiker);
                }
                //std::cout << "Synapse " << std::to_string(prePop) << " to " << std::to_string(post_population) << " at time " << std::to_string(this->infoGlobal->timeStep) << "\n";

                // for(signed int i_data = 0; i_data<currentBin.synapticState[post_population][prePop].size(); i_data++){
                // //std::cout << std::to_string(currentBin.synapticState[post_population][prePop][i_data]) << " ";
                // //std::cout << std::to_string(currentBin.synapticState[post_population][prePop][i_data]/currentBin.noRecordedSynapses[post_population][prePop]) << "\n";
                // }
            }
		}
	}
	RecordSynapseStates();
	RecordHeatmap();
	RecordRasterplot();
	RecordCurrents(synaptic_dV);
	//Record_Histogram(synaptic_dV);
	RecordAverages();
	RecordPotential();
	RecordCurrentContributions(synaptic_dV);
    if (this->heteroRecordPerSteps !=0 && (this->infoGlobal->timeStep % this->heteroRecordPerSteps) == 0) {
        //This condition will trigger if there is a "extra" number in the Parameters.txt and the step count is divisible by that number
        //1 for every step recording, 2 for every two steps, etc...
        RecordHeteroSynapses();
    	RecordHeteroSynapsesOverall();
    }
}

void Recorder::CloseStreams() {
    fileStreams.averagesFileStream.close();
    if (recordedHeatmap!=0){
        for (std::ofstream& individualFile : fileStreams.heatmapStreamVector) {
            individualFile.close();
        }
    } if (recordRasterPlot){
        fileStreams.rasterplotFileStream.close();
    } if (recordNeuronPotentials){
        fileStreams.potentialFileStream.close();
        fileStreams.currentsFileStream.close();
    } if (recordCurrentContributions){
        fileStreams.cCurrentsFileStream.close();
    } if (trackSynapses){
        fileStreams.synStatesFileStream.close();
    } if (recordHeteroSynapses){
        fileStreams.heteroSynFileStream.close();
        fileStreams.hSOverallFileStream.close();
    // } if (streamingNOutputBool){
    //     for (std::ofstream& stream : fileStreams.neuronOuputFileStreams){
    //         stream.close();
    // }
    }
}

void Recorder::WriteFinalDataFile(std::chrono::seconds setupTime, std::chrono::seconds simulationTime) {
    // std::chrono::time_point localTime = std::chrono::floor<std::chrono::seconds>(std::chrono::system_clock::now());//Does not work in GCC
    time_t timeInType = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

    std::stringstream outputString;
    // tm timeStruct = threadsafe::localtime(&timeInType);
    threadsafe::put_time(timeInType, "%d-%m-%Y %H:%M:%S", outputString);

    std::ofstream wParameterStream(GetParametersFilename(),std::ofstream::out | std::ofstream::app);
    wParameterStream <<  "#*****************************************************************\n";
    wParameterStream <<  "#Comp. finalized: " << outputString.rdbuf()<<"\tdd-mm-yyyy hh:mm:ss\n";
    wParameterStream <<  "#Set-up time:     " << setupTime.count() << " secs.\n";
    wParameterStream <<  "#Simulation time: " << simulationTime.count() << " secs.\n";
    // wParameterStream <<  "#This simulation ran on commit hash: "<<std::string(GIT_COMMIT)<<"\n";
    wParameterStream.close();
    std::cout << "Results saved in " << directoryPath << "\n";

}