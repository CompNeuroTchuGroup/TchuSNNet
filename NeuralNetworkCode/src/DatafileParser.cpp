//
// Created by Antoni Bertolin on 14.06.23
//
#include "DatafileParser.hpp"

DatafileParser::DatafileParser(Recorder &recorder) { // revised

  recorder.CloseStreams();
  parsingEnabled = recorder.parserEnabled;

  directoryPath   = recorder.directoryPath;
  simulationTitle = recorder.simulationTitle;
  for (PopInt neuronPop : std::ranges::views::iota(0, static_cast<PopInt>(recorder.noRasterPlotNeurons.size()))) {
    if (recorder.noRasterPlotNeurons.at(neuronPop) != 0) {
      neuronPopRasterIds.push_back(neuronPop);
    }
  }
  if (recorder.recordRasterPlot) {
    for (PopInt neuronPop : neuronPopRasterIds) {
      rasterPlotMetadata.push_back(recNeuronPopData(recorder.noRasterPlotNeurons.at(neuronPop), recorder.infoGlobal->dtTimestep,
                                                    static_cast<int>(recorder.infoGlobal->timeStep), neuronPop, recorder.infoGlobal->simulationTime));
    }
    fileNamesToParse.push_back(recorder.GetRasterplotFilename());
    fileTypes.push_back(RasterPlot);
  }
  for (std::string fileName : fileNamesToParse) {
    filesToParse.push_back(std::ifstream(fileName, std::ifstream::in));
  }
  // Here we have to access the recorder to obtain the filepaths/filenames and paths in general (and the metadata)
  // Then we have to convert the filepaths into the proper ifstreams
  // After this constructor it is no longer necessary to keep the NN object in scope
}

std::vector<FileEntry> DatafileParser::parseFileToEntries(std::ifstream &fileStream) { // revised
  // Here I have to parse each line (ignoring the first one with the column titles) and create FileEntries for each line
  std::vector<FileEntry> parsedFileEntries;
  std::string            entry;
  while (std::getline(fileStream, entry, '\n')) {
    if (entry.c_str()[0] == '#') {
      continue;
    } else if (entry.c_str()[0] == 'S') {
      continue;
    } else {
      parsedFileEntries.push_back(SplitStringToEntry(entry));
    }
  }
  fileStream.close();
  return parsedFileEntries;
}

void DatafileParser::setUpSpikeTimesVector(std::vector<std::vector<std::pair<std::vector<double>, std::pair<PopInt, NeuronInt>>>> &spikeTimesVector) {
  // population level
  for (PopInt recordedNeuronPop : std::ranges::views::iota(
           0, static_cast<PopInt>(
                  neuronPopRasterIds.size()))) { // To do proper indexing later, I will need a index-tracing function to get the index in the
                                                 // neuronPopRasterIds vector to index the larger vector properly when reading entries
    // neuron level
    std::vector<std::pair<std::vector<double>, std::pair<PopInt, NeuronInt>>> tempVector;
    for (NeuronInt neuron : std::ranges::views::iota(0, rasterPlotMetadata.at(recordedNeuronPop).noNeurons)) {
      // std::pair<std::vector<double>, std::pair<int, int>> tempElement;
      // tempElement.second=std::pair<int, int>(neuronPop, i);
      tempVector.push_back(std::pair<std::vector<double>, std::pair<PopInt, NeuronInt>>(
          {}, std::pair<PopInt, NeuronInt>(neuronPopRasterIds.at(recordedNeuronPop), neuron)));
    }
    spikeTimesVector.push_back(tempVector);
  }
}

std::vector<std::vector<std::pair<std::vector<double>, std::pair<PopInt, NeuronInt>>>> DatafileParser::parseSpikesToSpikeTimes(
    std::ifstream &fileStream) { // changing
  // Wrapper 1, assume the file does exist
  std::vector<std::vector<std::pair<std::vector<double>, std::pair<PopInt, NeuronInt>>>>
      spikeTimesVector; // vector of vectors of size equal to the number of neurons
  setUpSpikeTimesVector(spikeTimesVector);
  std::vector<FileEntry> parsedEntries{parseFileToEntries(fileStream)};
  if (parsedEntries.at(0).parameterValues.size() > 2) {
    throw "Error: the file is not properly formatted";
  }
  for (FileEntry &parsedEntry : parsedEntries) {
    int neuronPopIndex = static_cast<PopInt>(std::distance(
        neuronPopRasterIds.begin(), std::find(neuronPopRasterIds.begin(), neuronPopRasterIds.end(),
                                              std::stoi(parsedEntry.parameterValues.at(1))))); // How do we do exception management here?
    if (neuronPopIndex == static_cast<PopInt>(std::distance(neuronPopRasterIds.begin(), neuronPopRasterIds.end()))) {
      throw "Indexing error: the neuron population was not found";
    }
    spikeTimesVector.at(neuronPopIndex).at(std::stoi(parsedEntry.parameterValues.at(0))).first.push_back(std::stod(parsedEntry.parameterName));
  }

  // Here we can use the FileEntry struct functions to parse a line, then iterate over parameterValues(two different loops). Name will be the
  // timestep, luckily enough. If 1 is found, timestep added to the neuron's vector. It will be necessary to properly use the indexes to achieve this
  // (obtained in the  function above) Use this return to call the write function
  return spikeTimesVector;
}

void DatafileParser::writeSpikeTimesFile(std::vector<std::vector<std::pair<std::vector<double>, std::pair<PopInt, NeuronInt>>>> parsedData,
                                         std::string wfilePath, std::vector<recNeuronPopData> metadataVec) { // later {
  // Wrapper 2
  std::ofstream stream(wfilePath);
  // Metadata loop
  for (recNeuronPopData singleMetadata : metadataVec) {
    stream << "M_" << std::to_string(singleMetadata.neuronPopId) << '=' << std::to_string(singleMetadata.noNeurons) << ','
           << std::to_string(singleMetadata.dtTimestep) << ',' << std::to_string(singleMetadata.totalTimesteps) << ','
           << std::to_string(singleMetadata.neuronPopId) << ',' << std::to_string(singleMetadata.simTime) << '\n';
  }
  // Population loop
  for (std::vector<std::pair<std::vector<double>, std::pair<PopInt, NeuronInt>>> &population : parsedData) {
    // neuron loop
    for (std::pair<std::vector<double>, std::pair<PopInt, NeuronInt>> &neuron : population) {
      stream << "N_" << std::to_string(neuron.second.first) << '_' << std::to_string(neuron.second.second) << "=";
      for (double spiketime : neuron.first) {
        stream << std::to_string(spiketime) << ",";
      }
      stream << "\n";
    }
  }
  stream.close();
  // Here I will have to first create the ofstream. Then basically iterate over the vector and print each vector with a comma separator and write in
  // desired format: M(etadata)=noNeurons,dt,totaltimesteps,neuronPopId N(euronpop)_1=spiketime1,spiketime2,
}

void DatafileParser::parse() {
  // Parent wrapper of wrappers 1 and 2
  if (parsingEnabled) {
    int index{};
    for (std::ifstream &fileStream : filesToParse) {
      struct stat buffer;
      if (stat(fileNamesToParse.at(index).c_str(), &buffer) != 0) {
        std::cout << "*************************\n";
        std::cout << fileNamesToParse.at(index) << " file does not exist\n";
        std::cout << "*************************\n";
        throw "File does not exist";
      }
      if (neuronPopRasterIds.size() != 0 && (fileTypes.at(index) == RasterPlot)) {
        // Now we assume the file exist
        writeSpikeTimesFile(parseSpikesToSpikeTimes(fileStream), this->getParsedSpikeTimesFilePath(), rasterPlotMetadata);
      } else {
        throw "Logical error in DatafileParser.cpp";
      }
      index++;
    }
  }
  std::cout << "\nParsing operations are finished.\n";
}

recNeuronPopData::recNeuronPopData(NeuronInt noNeurons, double dtTimestep, int totalTimesteps, PopInt neuronPopId, double simTime)
    : noNeurons{noNeurons}, dtTimestep{dtTimestep}, simTime{simTime}, totalTimesteps{totalTimesteps}, neuronPopId{neuronPopId} {
}
