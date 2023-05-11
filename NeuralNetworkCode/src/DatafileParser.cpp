#include "DatafileParser.hpp"

DatafileParser::DatafileParser(AdvancedRecorder& recorder) //revised
{
    recorder.CloseStreams();
    parsingEnabled=recorder.parserEnabled;

    directoryPath=recorder.directoryPath;
    title=recorder.title;
    for (int i = 0; i < recorder.noRasterPlotNeurons.size(); i++){
        if (recorder.noRasterPlotNeurons[i] != 0){
            neuronPopRasterIds.push_back(i);
        }
    }
    if(recorder.noRasterPlotNeurons.sum() != 0){
        for (int neuronPop : neuronPopRasterIds){
            rasterPlotMetadata.emplace_back(recNeuronPopData(recorder.noRasterPlotNeurons[neuronPop],recorder.info->dt,static_cast<int>(recorder.info->time_step), neuronPop, recorder.info->simulationTime));
        }
        fileNamesToParse.emplace_back(recorder.GetRasterplotFilename());
        fileTypes.push_back(RasterPlot);
    }
    for (std::string fileName : fileNamesToParse){
            filesToParse.emplace_back(std::ifstream(fileName, std::ifstream::in));
    }
    //Here we have to access the recorder to obtain the filepaths/filenames and paths in general (and the metadata)
    //Then we have to convert the filepaths into the proper ifstreams
    //After this constructor it is no longer necessary to keep the NN object in scope
}

std::vector<FileEntry> DatafileParser::parseFileToEntries(std::ifstream& fileStream) //revised
{    
    //Here I have to parse each line (ignoring the first one with the column titles) and create FileEntries for each line
    std::vector<FileEntry> parsedFileEntries;
    for (std::string entry; std::getline(fileStream, entry);){
        if(entry.c_str()[0] == '#'){
            continue;
        } else if (entry.c_str()[0] == 'S'){
            continue;
        } else {
            parsedFileEntries.emplace_back(stringToFileEntry(std::move(static_cast<std::string>(entry))));
        }
    }
    fileStream.close();
    return parsedFileEntries;
}

void DatafileParser::setUpSpikeTimesVector(std::vector<std::vector<std::pair<std::vector<double>, std::pair<int, int>>>> &spikeTimesVector)
{
    //population level
    for (int index = 0; index<neuronPopRasterIds.size(); index++){ //To do proper indexing later, I will need a index-tracing function to get the index in the neuronPopRasterIds vector to index the larger vector properly when reading entries
    //neuron level
        std::vector<std::pair<std::vector<double>, std::pair<int, int>>> tempVector;
        for (int i = 0; i<rasterPlotMetadata.at(index).noNeurons; i++){
            // std::pair<std::vector<double>, std::pair<int, int>> tempElement;
            // tempElement.second=std::pair<int, int>(neuronPop, i);
            tempVector.push_back(std::pair<std::vector<double>, std::pair<int, int>>({}, std::pair<int, int>(neuronPopRasterIds.at(index), i)));
        }
        spikeTimesVector.push_back(tempVector);
    }
}

std::vector<std::vector<std::pair<std::vector<double>, std::pair<int, int>>>> DatafileParser::parseSpikesToSpikeTimes(std::ifstream &fileStream) // changing
{
    //Wrapper 1, assume the file does exist
    std::vector<std::vector<std::pair<std::vector<double>, std::pair<int, int>>>> spikeTimesVector; //vector of vectors of size equal to the number of neurons
    setUpSpikeTimesVector(spikeTimesVector);
    std::vector<FileEntry> parsedEntries{parseFileToEntries(fileStream)};
    if (parsedEntries[0].values.size() > 2) {
        std::cout << "Error: the file is not properly formatted"<<"\n";
        throw;
    }
    for (FileEntry& parsedEntry : parsedEntries){
        int neuronPopIndex = static_cast<int>(std::distance(neuronPopRasterIds.begin(), std::find(neuronPopRasterIds.begin(),neuronPopRasterIds.end(), std::stoi(parsedEntry.values.at(1))))); //How do we do exception management here?
        if (neuronPopIndex == static_cast<int>(std::distance(neuronPopRasterIds.begin(), neuronPopRasterIds.end()))){
            throw;
            std::cout<<"Indexing error: the neuron population was not found"<<"\n";
        }
        spikeTimesVector.at(neuronPopIndex).at(std::stoi(parsedEntry.values.at(0))).first.push_back(std::stod(parsedEntry.name));
    }

    //Here we can use the FileEntry struct functions to parse a line, then iterate over values(two different loops). Name will be the timestep, luckily enough. 
    //If 1 is found, timestep added to the neuron's vector.
    //It will be necessary to properly use the indexes to achieve this (obtained in the  function above)
    //Use this return to call the write function
    return spikeTimesVector;
}



void DatafileParser::writeSpikeTimesFile(std::vector<std::vector<std::pair<std::vector<double>, std::pair<int, int>>>> parsedData, std::string wfilePath, std::vector<recNeuronPopData> metadataVec)//later
{
    //Wrapper 2
    std::ofstream stream(wfilePath);
    //Metadata loop
    for (recNeuronPopData singleMetadata : metadataVec){
        stream<<"M_"<<std::to_string(singleMetadata.neuronPopId)<<'='<<std::to_string(singleMetadata.noNeurons)<<','<<std::to_string(singleMetadata.dt)<<','<<std::to_string(singleMetadata.totalTimesteps)<<','<<std::to_string(singleMetadata.neuronPopId)<<','<<std::to_string(singleMetadata.simTime)<<'\n';
    }
    //Population loop
    for (std::vector<std::pair<std::vector<double>, std::pair<int, int>>>& population : parsedData){
        //neuron loop
        for (std::pair<std::vector<double>, std::pair<int, int>>& neuron : population){
            stream<<"N_"<<std::to_string(neuron.second.first)<<'_'<<std::to_string(neuron.second.second)<<"=";
            for (double spiketime : neuron.first){
                stream<<std::to_string(spiketime)<<",";
            }
            stream<<"\n";
        }
    }
    stream.close();
    //Here I will have to first create the ofstream. Then basically iterate over the vector and print each vector with a comma separator and write in desired format:
    //M(etadata)=noNeurons,dt,totaltimesteps,neuronPopId
    //N(euronpop)_1=spiketime1,spiketime2,
}

void DatafileParser::parse()
{
    //Parent wrapper of wrappers 1 and 2
    if (parsingEnabled){
        int index{};
        for (std::ifstream& fileStream : filesToParse){
            struct stat buffer;
            if(stat(fileNamesToParse.at(index).c_str(),&buffer) != 0){
                std::cout << "*************************\n";
                std::cout <<fileNamesToParse.at(index)<<" file does not exist\n";
                std::cout << "*************************\n";
                throw;
            }
            if(neuronPopRasterIds.size() !=0 && (fileTypes.at(index)==RasterPlot)){
                //Now we assume the file exist
                writeSpikeTimesFile(parseSpikesToSpikeTimes(fileStream), this->GetParsedSpikeTimesFilePath(), rasterPlotMetadata);
            } else {
                std::cout<<"\nSomething went wrong\n";
                throw;
            }
            index++;
        }
    }
    std::cout<<"\nParsing operations are finished.\n";
}

recNeuronPopData::recNeuronPopData(int noNeurons, double dt, int totalTimesteps,int neuronPopId, double simTime): noNeurons{noNeurons}, dt{dt}, totalTimesteps{totalTimesteps}, neuronPopId{neuronPopId}, simTime{simTime}
{
}
