#ifndef Recorder_HPP
#define Recorder_HPP

// #define _CRT_SECURE_NO_WARNINGS
#include "./GlobalFunctions.hpp"
#include "NeuronPopSample.hpp"
#include "Stimulus/Stimulus.hpp"
#include "SynapseSample.hpp"
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#ifdef _WIN32
#include <direct.h>
// #elif __APPLE__
#endif

#include <cstring>
#include <fstream>
#include <iostream>

#include <chrono>

struct binData {
  std::vector<double> potential;   // used to compute average potential of neurons [per population]
  std::vector<double> spikerRatio; // number of spiked neurons [per population]

  std::vector<double>              externalCurrent;  // average external input current [to population totalNeuronPopulation] //to vector
  std::vector<std::vector<double>> synapticCurrents; // average current [to population P1][from population P2] //to vector
  // std::valarray<double>  totalCurrent_mean;             // mean of the total input current     [to population totalNeuronPopulation]

  std::vector<std::vector<double>> neuronTotalCurrentMean;   // mean of the total input current     [to neuron N]
  std::vector<double>              totalCurrentSquared_mean; // mean of the squared total input current to each neuron [averaged over population Ps]
  // std::valarray<int>               spiker_pop; //Member is currently not in use

  /** synaptic statistics per time step for every postsynaptic and every
   * presynaptic population and every synapse specific data column.
   */
  std::vector<std::vector<std::vector<double>>> synapticState; // synaptic data of synapses [to population P1][from population P2][data entry j] //sum
                                                               // of vectors, valarray is valid (in last level)
  std::vector<std::vector<double>>      heatmap;               // firing rate [of population i][in each pixel]
  std::vector<std::vector<signed long>> noRecordedSynapses;
};

class Recorder {
  friend class DatafileParser;

protected:
  static constexpr size_t recordBufferSize{131072}; // Do not change, optimal buffer size for data recording

  GlobalSimInfo *infoGlobal;

  std::string simulationTitle;
  std::string nonIterateTitle;
  std::string directoryPath;
  std::string dateTime;

  std::shared_ptr<const NeuronPopSample> neurons;
  std::shared_ptr<const SynapseSample>   synapses;
  std::shared_ptr<const Stimulus>        stimulus;

  int       timeStepsPerBin{};
  NeuronInt noNeuronsConnectivity{};
  NeuronInt noNeuronsDelay{};
  NeuronInt noNeuronsJPot{};
  // AdvRecorder parameters
  bool parserEnabled{false};
  bool trackSynapses{false}; //, writeHistogram;
  int  recordedHeatmap{0};
  int  initialCurrent{0};
  int  startRecordingTime{0};

  std::vector<NeuronInt>           neuronPotentialsToRecord;
  bool                             recordNeuronPotentials{false};
  std::vector<NeuronInt>           noRasterPlotNeurons;
  bool                             recordRasterPlot{false};
  std::vector<NeuronInt>           currentContributionsToRecord;
  bool                             recordCurrentContributions{false};
  std::vector<std::vector<double>> densimap; // Number of neurons [of population i][in each pixel] //Is densimap necessary?
  std::vector<double>              currentContrBin;
  // statistics per time step
  binData currentBin;

  std::vector<std::pair<NeuronInt, signed long>> heteroSynTracker; //
  bool                                           recordHeteroSynapses{false};
  RecorderOpenStreams                            fileStreams;

  TStepInt heteroRecordPerSteps{10};

  void ResetStatistics(); // Resets all containers.

  // Record functions
  void RecordPotential();
  void RecordRasterplot();
  void RecordCurrents(const std::vector<std::vector<double>> &synaptic_dV);
  void RecordAverages();
  void RecordSynapseStates();
  void RecordHeatmap();
  void RecordCurrentContributions(const std::vector<std::vector<double>> &synaptic_dV);
  void RecordHeteroSynapses();
  void RecordHeteroSynapsesOverall();

  // WriteDataHeader functions
  void WriteDataHeaderCurrents();
  void WriteDataHeaderRasterplot();
  void WriteDataHeaderAverages();
  void WriteDataHeaderSynapseStates();
  void WriteDataHeaderPotential();
  void WriteDataHeaderHeatmap();
  void WriteDataHeaderCurrentsContribution();
  void WriteDataHeaderHeteroSynapses();
  void WriteDataHeaderHeteroSynapsesOverall();

  void SetNoRasterplotNeurons(const std::vector<std::string> &parameterValues);
  void SetNoRecordedNeuronPotentials(const std::vector<std::string> &parameterValues);
  void SetNoCurrentContribution(const std::vector<std::string> &parameterValues);
  void SetNoRecordedHeteroSynapticProfilesPerTrackedNeuronPerPop(const std::vector<std::string> &parameterValues);

  void SetAveragingSteps(double secondsPerBin);
  void BindNoHeteroSynapsesPerPop(PopInt neuronPop);

  void AllocateAndAssignStreamBuffer(std::ofstream &outputStream);

public:
  Recorder(const std::shared_ptr<NeuronPopSample> neurons, const std::shared_ptr<SynapseSample> synapses, const std::shared_ptr<Stimulus> stimulus,
           std::string baseDirectory, std::vector<FileEntry> inputParameters, std::string titleString, std::string nonIterateTitle,
           GlobalSimInfo *infoGlobal);
  ~Recorder() = default;

  void WriteDataHeader();
  void WriteFinalDataFile(std::chrono::seconds setupTime, std::chrono::seconds simulationTime);
  void Record(const std::vector<std::vector<double>> &synapticInput);

  //***** Set-Functions *****

  void SetFilenameDate();

  //**** Get-Functions *****
  int         GetAveragingSteps() const { return timeStepsPerBin; }
  std::string GetDirectoryPath() const { return this->directoryPath + simulationTitle; }
  std::string GetDataFilename() const { return this->directoryPath + simulationTitle + "_Data.dat"; }
  std::string GetParametersFilename() const { return this->directoryPath + simulationTitle + "_Parameters.txt"; }
  std::string GetParameterOptionsFilename() const { return this->directoryPath + "ParameterOptions.txt"; }
  std::string GetConnectivityFilename() const { return this->directoryPath + simulationTitle + "_Connectivity_Matrix"; }
  std::string GetDelayFilename() const { return this->directoryPath + simulationTitle + "_DelayConnectivity_Matrix"; }
  std::string GetJPotFilename() const { return this->directoryPath + simulationTitle + "_JPotConnectivity_Matrix"; }
  std::string GetHeatmapFilename() const { return this->directoryPath + simulationTitle + "_HeatmapRate_Pop"; }
  std::string GetTitle() const { return simulationTitle; }
  // Legacy advanced recorder: get functions
  std::string GetRasterplotFilename() const { return this->directoryPath + simulationTitle + "_Rasterplot.dat"; }
  std::string GetCurrentsFilename() const { return this->directoryPath + simulationTitle + "_Currents.dat"; }
  std::string GetPotentialFilename() const { return this->directoryPath + simulationTitle + "_Potential.dat"; }
  std::string GetSynapseStateFilename() const { return this->directoryPath + simulationTitle + "_Synapses.dat"; }
  std::string GetCurrentCrontributionFilename() const { return this->directoryPath + simulationTitle + "_CurrentContribution.dat"; }
  std::string GetHeteroSynapseStateFilename() const { return this->directoryPath + simulationTitle + "_HeteroSynapses.dat"; }
  std::string GetOverallHeteroSynapseStateFilename() const { return this->directoryPath + simulationTitle + "_OverallHS.dat"; }
  // std::string GetHeteroBranchedSynapseStateFilename() const { return this->directoryPath + simulationTitle + "_BranchedHS.dat"; }
  std::string GetSteadyStatesFilename() const { return this->directoryPath + simulationTitle + "_steadyStates.dat"; }

  void LoadParameters(const std::vector<FileEntry> &input);
  void SaveParameters(std::ofstream &stream) const;

  void CloseStreams();

  void WriteHeader(std::ofstream &fileStream) const;
  void WriteConnectivity() const;
  void WriteDistributionD() const;
  void WriteDistributionJ() const;
  void WriteSteadyStates() const;

  void MakeInputCopies(const std::string &filename) const;
};

#endif /* Recorder_HPP */
