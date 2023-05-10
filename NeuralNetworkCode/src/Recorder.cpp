
#include "Recorder.hpp"

void Recorder::SetFilenameDate(){
	bool Windows = false;

    #if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
        Windows = true;
    #endif
	
    time_t      t = time(NULL);
    tm*         timePtr = localtime(&t);

    dateTime =  std::to_string(timePtr->tm_year + 1900) + "_" +
    std::to_string(timePtr->tm_mon+1) + "_" +
    std::to_string(timePtr->tm_mday) + "_" +
    std::to_string(timePtr->tm_hour) + "-" +
    std::to_string(timePtr->tm_min) + "-" +
    std::to_string(timePtr->tm_sec);
	
    if (Windows)
        directoryPath += title + "_" + dateTime + "\\";
    else
        directoryPath += title + "_" + dateTime +  "/";

#ifdef _WIN32
    _mkdir(directoryPath.c_str()); //0744
#elif __APPLE__
    mkdir(directoryPath.c_str(),0744);
#elif __linux__
    mkdir(directoryPath.c_str(),0744);
#endif


}

Recorder::Recorder(NeuronPopSample *ns,SynapseSample *syn,Stimulus *sti,std::string baseDir,std::vector<std::string> *input,std::string str_t,GlobalSimInfo * info){
    this->info            = info;
    neurons               = ns;
    synapses              = syn;
    stimulus              = sti;
    this->directoryPath   = baseDir;
    title                 = str_t;
    noNeuronsConnectivity = 0;
    noNeuronsDelay        = 0;
    noNeuronsJPot         = 0;
}

void Recorder::WriteConnectivity(){
    this->synapses->WriteConnectivity(GetConnectivityFilename(),noNeuronsConnectivity);
}

void Recorder::WriteDistributionD(){
    this->synapses->WriteDistributionD(GetDelayFilename(),noNeuronsDelay);
}

void Recorder::WriteDistributionJ(){
    this->synapses->WriteDistributionJ(GetJPotFilename(),noNeuronsJPot);
}

void Recorder::SaveParameters(std::ofstream * stream){
    *stream <<  "#*************************************************\n";
    *stream <<  "#************** Recorder Parameters **************\n";
    *stream <<  "#*************************************************\n";
    *stream <<  "recorder_type\t\t\t\t\t\t" << GetType() <<  "\n";
    *stream <<  "recorder_noNeuronsConnectivity\t\t" << std::to_string(noNeuronsConnectivity) << "\t\t\t\t\t#saves connectivity matrices for the first x neurons of each populations\n";
    *stream <<  "recorder_noNeuronsDelay\t\t\t\t" << std::to_string(noNeuronsDelay) << "\t\t\t\t\t#saves delay connectivity matrices for the first x neurons of each populations\n";
    *stream <<  "recorder_noNeuronsJPot\t\t\t\t" << std::to_string(noNeuronsJPot) << "\t\t\t\t\t#saves Jpot connectivity matrices for the first x neurons of each populations\n";
}


void Recorder::LoadParameters(std::vector<std::string> *input){
    std::string              name,token;
    std::vector<std::string> values;

    for(std::vector<std::string>::iterator it = (*input).begin(); it != (*input).end(); ++it) {

        SplitString(&(*it),&name,&values);

        if(name.find("recorder_noNeuronsConnectivity") != std::string::npos)
            noNeuronsConnectivity = std::stoi(values.at(0));
        else if(name.find("recorder_noNeuronsDelay") != std::string::npos)
            noNeuronsDelay = std::stoi(values.at(0));
        else if(name.find("recorder_noNeuronsJPot") != std::string::npos)
            noNeuronsJPot = std::stoi(values.at(0));

    }
}


void Recorder::SetAveraging_steps(double rec_dt) {
    double dt = info->dt;
    if(rec_dt < dt)
        averaging_steps = 1;
    else
        averaging_steps = static_cast<int>(std::round(rec_dt / dt));
}

void Recorder::WriteHeader(std::ofstream * stream){
    *stream <<  "#*****************************************************************\n";
    *stream <<  "# Time and Title: " << dateTime << " -- " << title << "\n";
    *stream <<  "#*****************************************************************\n";
}

void Recorder::makeInputCopy(const std::string& inputFile) {
    std::ifstream  src(inputFile, std::ios::binary);
    std::ofstream  dst(this->directoryPath + "input_copy_parameters.txt",   std::ios::binary);
    dst << src.rdbuf();
}
