#ifndef _DATA_FILE_PARSER_HPP_
#define _DATA_FILE_PARSER_HPP_
#define _MAX_CHARACTERS_PER_LINE 30002//INT_MAX/2 // This may require adjustment

#include "AdvancedRecorder.hpp"
#include "GlobalFunctions.hpp"

#include <stdlib.h>
#include <errno.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

struct recNeuronPopData{
    int noNeurons{};
    double dt{};
    double simTime{};
    int totalTimesteps{};
    int neuronPopId{};
    recNeuronPopData(int noNeurons, double dt, int totalTimesteps,int neuronPopId, double simTime);
};
enum fileType{
    RasterPlot, GeneralData, Currents, SynapseStates, Potential, Heatmap, CurrentsContribution, HeteroSynapseStates, HeteroSynapseOverall
};

class DatafileParser {

protected:

    bool parsingEnabled{false};
    //Bools will be necessary to check which parsing is needed, or whether it is needed or not
    //Add here more if there are more files to be parsed at the end of the simulation

    std::string directoryPath;
    std::string title;
    std::string extension {".dat"};

    std::vector<int> neuronPopRasterIds;

    std::vector<std::string> fileNamesToParse; //Here we store the file names to open. A similar path will be used for writing, adding _parsed to the filename
    std::vector<std::ifstream> filesToParse; //Here we store the streams to read. Streams to write can be  opened and closed locally.
    std::vector<fileType> fileTypes; // To have a way to know the filetype we are parsing

    std::vector<recNeuronPopData> rasterPlotMetadata; //Here we s tore any extra information to write in the file that may be necessary for reading the datafile in the future.

public:
    //Constructor/destructor
    explicit DatafileParser(AdvancedRecorder& recorder);
    ~DatafileParser()=default;

    std::string GetParsedSpikeTimesFilePath() {return directoryPath + title + "_ParsedSpikeTimes" + extension;}

    void setUpSpikeTimesVector(std::vector<std::vector<std::pair<std::vector<double>, std::pair<int, int>>>>& spikeTimesVector);
    std::vector<std::vector<std::pair<std::vector<double>, std::pair<int, int>>>> parseSpikesToSpikeTimes(std::ifstream& fileStream); //Create a vector of vectors (one per neuron) of spike times from the read file

    std::vector<FileEntry> parseFileToEntries(std::ifstream& fileStream);

    void writeSpikeTimesFile(std::vector<std::vector<std::pair<std::vector<double>, std::pair<int, int>>>> parsedData, std::string wfilePath, std::vector<recNeuronPopData> metadataVEC); //
    //bool closeOpenFile(int index);

    void parse();

};



#endif