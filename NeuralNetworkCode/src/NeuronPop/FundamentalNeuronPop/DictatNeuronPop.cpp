//
// Created by Antoni Bertolin on 14.06.23
//
#include "./DictatNeuronPop.hpp"
#include "DictatNeuronPop.hpp"

DictatNeuronPop::DictatNeuronPop(GlobalSimInfo* infoGlobal, NeuronInt neuronID) : NeuronPop (infoGlobal, neuronID) {
    //Constructor
    inputFileAddress = infoGlobal->pathToInputFile + "DictatNeuronPop_"+std::to_string(neuronID)+"_spikers.txt";
    inputStreamFile.open(inputFileAddress, std::ifstream::in);
    //Then just open and read, with bools for raw spikers or instruction-based
}

void DictatNeuronPop::LoadParameters(const std::vector<FileEntry>& neuronParameters)
{
    NeuronPop::LoadParameters(neuronParameters);//For now this is commented out to avoid confusion, and prevent from recording output in DictatNeuronPop

    for(auto& [parameterName, parameterValues] : neuronParameters) {

        if(parameterName.find("inputFile") != std::string::npos){
            if(parameterValues.at(0).find("spiker") != std::string::npos){
                spikerListFiringBool=true;
            } else if (parameterValues.at(0).find("instruction") != std::string::npos){
                instructionFiringBool=true;
                inputInstructions.resize(noNeurons);
                activeInstructions.resize(noNeurons);
                ReadInstructionsFromFile();
            } else {
                throw "Non-existant instruction type"; //A non-existant type should throw
            }
        }else if(parameterName.find("poissonLike") != std::string::npos){
            if(parameterValues.at(0).find("true") != std::string::npos){
                poissonLikeFiringBool=true;
            } else if (parameterValues.at(0).find("false") != std::string::npos){
                poissonLikeFiringBool=false;
            } else {
                throw "Non-existant bool in DictatNeuronPop"; //A non-existant type should throw
            }
        }  
    }
}

void DictatNeuronPop::SaveParameters(std::ofstream& wParameterStream) const {

    std::string idString = "neurons_" + std::to_string(GetId());

    NeuronPop::SaveParameters(wParameterStream);

    wParameterStream <<  idString + "_inputFile\t\t\t";
    if (spikerListFiringBool){
        wParameterStream<<"spiker";
    } else if (instructionFiringBool){
        wParameterStream<<"instruction";
    } else {
        wParameterStream<<"None";
    }
    wParameterStream <<  "\t#Write 'instruction' for instruction-based, and 'spiker' for spiker-based \n";

    wParameterStream <<  idString + "_poissonLike\t\t\t";
    if (poissonLikeFiringBool){
        wParameterStream<<"true";
    } else {
        wParameterStream<<"false";
    }
    wParameterStream <<  "\t#Write true vs false to indicate if the firing reproduced from instructions is poisson-like or periodic. \n";
}

void DictatNeuronPop::Advect(const std::vector<double>& synaptic_dV) {
    spikerNeuronsPrevdt = spikerNeurons;
    spikerNeurons.clear();
    if (spikerListFiringBool){
        ReadSpikersFromFile();
    } else if(instructionFiringBool){
        if (poissonLikeFiringBool){
            GeneratePoissonSpikersFromInstructions();
        } else {
            GenerateRegularSpikersFromInstructions();
        }
    } else {
        throw "Logical error in DictatNeuronPop";
    }

}

void DictatNeuronPop::ReadInstructionsFromFile() {
    /*Instruction files for a population should be written in the following format:
    > neuronid starttime1 endtime1 fire_every_n_steps
    And ALWAYS time ordered with correct intervals (wrong starttimes will essentially change the first AP firing timestep or total phase)->(timeste-starttime)%frequency
    */
   //Runs only once in LP
   std::string  line;
    while(std::getline(inputStreamFile,line, '\n')){
        if(line[0] == '#'){
            continue;
        }
        else if (line[0] == '>') {
            FileEntry entry{ SplitStringToEntry(line) };
            inputInstructions.at(std::stoi(entry.parameterValues.at(0))).emplace_back(Instruction(entry, infoGlobal->dtTimestep));
        } else if (line[0] == ' ' || line[0] == '\n') {
            continue;
        } else {
            std::cout<<"Reading error. The unexpected read input was: "<< line <<"\n";
            throw strncat(line.data(), " -> was unexpected reading input.Reading error.", 48);
        }
    }
    for (std::vector<Instruction>& neuronInstructions: inputInstructions){
        neuronInstructions.back().last=true;
    }
}

void DictatNeuronPop::GenerateRegularSpikersFromInstructions()
{
    // //loop option 1 (inefficient but less exception prone)
    // for (int neuronId = 0; neuronId<static_cast<int>(inputInstructions.size()); neuronId++){
    //     for (Instruction& instruction : inputInstructions.at(neuronId)){
    //         if ((!instruction.completed) && (instruction.startTimeStep>infoGlobal->timeStep)){ //Here not >= to avoid a spike of frequency when changing instructions
    //             if (((infoGlobal->timeStep-instruction.startTimeStep)%instruction.fireEveryNSteps)==0){
    //                 spiker.push_back(neuronId);
    //                 if (instruction.endTimeStep>=infoGlobal->timeStep){
    //                 instruction.completed=true;
    //             }
    //             break;
    //             } 
    //         }
    //     }
    // }
    //loop option 2 (efficient but probably worse in exception handling)
    for (NeuronInt neuronId : std::ranges::views::iota(0,noNeurons)){
        Instruction& instruction = inputInstructions.at(neuronId).at(activeInstructions.at(neuronId));
        if (instruction.endTimeStep<=infoGlobal->timeStep){
            if (!(instruction.last)){
                activeInstructions.at(neuronId)++;
            }
            instruction.completed=true;
        }
        if (!instruction.off && ((infoGlobal->timeStep-instruction.startTimeStep)%instruction.fireEveryNSteps)==0){
            //Instruction.off condition is there to avoid doing a modulus with zero
            if (instruction.completed){continue;}
            if (instruction.startTimeStep<infoGlobal->timeStep) {spikerNeurons.push_back(neuronId);}
        }
    }
}

void DictatNeuronPop::GeneratePoissonSpikersFromInstructions() {
        for (NeuronInt neuronId : std::ranges::views::iota(0,noNeurons)){
        Instruction& instruction = inputInstructions.at(neuronId).at(activeInstructions.at(neuronId));
        if (instruction.endTimeStep<=infoGlobal->timeStep){
            if (!(instruction.last)){
                activeInstructions.at(neuronId)++;
            }
            instruction.completed=true;
        }
        if (uniformDistribution(generator)<instruction.firingProbability){ //Here because there is no modulus, there is no need for checking the instruction.off condition
            if (instruction.completed){continue;}
            if (instruction.startTimeStep<infoGlobal->timeStep) {spikerNeurons.push_back(neuronId);}
        }
    }
}

void DictatNeuronPop::ReadSpikersFromFile()
{
    /*Spikers files for a population should be written in the following format:
    > time neuron1, neuron2, neuron 4.... (neurons that spike in timeStep)
    *NEVER PUT COMMENTS OR EXTRA_LINES IN THIS TYPE OF FILE!
    *And be careful with your dt, as this will alter behaviour
    *This file should be generated automatically and tested against Test11 (deorecated)
    *This function runs once per timestep
    */
    std::string line;
    while(std::getline(inputStreamFile, line, '\n')){
        if (line[0]=='>'){
           FileEntry s_entry {SplitStringToEntry(std::move(static_cast<std::string>(line)))};
            long timeStep{std::lround(std::stod(s_entry.parameterValues.at(0))/infoGlobal->dtTimestep)};
            if (timeStep != infoGlobal->timeStep){
                throw "There was a reading alignment error";
            } else {
                for (std::string value : s_entry.parameterValues){
                    spikerNeurons.push_back(std::stoi(value));
                }
            };
        } else {
            throw strncat(line.data(), " -> was unexpected reading input.Reading error.",48);
        }
    }
}

Instruction::Instruction(FileEntry inputEntry, double dtTimestep): neuronId{std::stoi(inputEntry.parameterValues.at(0))}, startTimeStep{std::lround(std::stod(inputEntry.parameterValues.at(1))/dtTimestep)}, endTimeStep{std::lround(std::stod(inputEntry.parameterValues.at(2))/dtTimestep)}, frequency{std::stod(inputEntry.parameterValues.at(3))}, firingProbability{frequency*dtTimestep} {
    //Constructor
    if (frequency < std::numeric_limits<double>::epsilon()){ //Zero comparison to avoid division by zero
        this->off=true;
    } else {
        fireEveryNSteps=std::lround((1/frequency)/dtTimestep); //Conversion from frequency to timestep period.
        if (fireEveryNSteps==0){ //If the frequency is close
            std::cout<<"\n"<<"EXCEPTION: YOU CHOSE A FREQUENCY THAT IS TOO HIGH FOR NEURON "<<std::to_string(neuronId)<<"\n\n\n"<<"**********************************";
            throw "EXCEPTION: YOU CHOSE A FREQUENCY THAT IS TOO HIGH FOR CURRENT DT IN DICTAT INPUT FILE";
        }
    }
}
