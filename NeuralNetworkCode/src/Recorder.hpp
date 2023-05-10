#ifndef Recorder_HPP
#define Recorder_HPP

#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cstring>
#include <fstream>
#include "NeuronPopSample.hpp"
#include "SynapseSample.hpp"
#include "Stimulus/Stimulus.hpp"
#include <chrono>
#include <ctime>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef _WIN32
    #include <direct.h>
//#elif __APPLE__
#endif

class Recorder
{
protected:

    GlobalSimInfo * info;

    std::string title;
    std::string directoryPath;
    std::string dateTime;

    NeuronPopSample     *neurons;
    SynapseSample       *synapses;
    Stimulus            *stimulus;

    int         averaging_steps;
    long        noNeuronsConnectivity;
    long        noNeuronsDelay;
    long        noNeuronsJPot;

public:
    Recorder(NeuronPopSample *ns,SynapseSample *syn,Stimulus *sti,std::string baseDir,std::vector<std::string> *input,std::string str_t,GlobalSimInfo * info);
    virtual ~Recorder()=default;

    virtual void WriteDataHeader() = 0;
	virtual void writeFinalDataFile(double comp_time) = 0;
	virtual void Record(std::vector<std::vector<double>> * synapticInput) = 0;

    //***** Set-Functions *****
    void SetAveraging_steps(double rec_dt);
    void SetFilenameDate();

    //**** Get-Functions *****
    int     GetAveragingSteps()   {return averaging_steps;}
    virtual std::string GetType() = 0;
    std::string GetDirectoryPath() {return this->directoryPath + title;}
    std::string GetDataFilename()        {return this->directoryPath + title + "_Data.dat";}
    std::string GetParametersFilename()  {return this->directoryPath + title + "_Parameters.txt";}
    std::string GetParameterOptionsFilename()  {return this->directoryPath + "ParameterOptions.txt";}
    std::string GetConnectivityFilename(){return this->directoryPath + title + "_Connectivity_Matrix";}
    std::string GetDelayFilename(){return this->directoryPath + title + "_DelayConnectivity_Matrix";}
    std::string GetJPotFilename(){return this->directoryPath + title + "_JPotConnectivity_Matrix";}
	std::string GetHeatmapFilename() { return this->directoryPath + title + "_HeatmapRate_Pop"; }
    std::string GetTitle()               {return title;}

    virtual void LoadParameters(std::vector<std::string> *input) = 0;
    virtual void SaveParameters(std::ofstream * stream);

    virtual void CloseStreams() = 0;

    void WriteHeader(std::ofstream * stream);
    void WriteConnectivity();
    void WriteDistributionD();
    void WriteDistributionJ();

    void makeInputCopy(const std::string&);

};

#endif /* Recorder_HPP */
