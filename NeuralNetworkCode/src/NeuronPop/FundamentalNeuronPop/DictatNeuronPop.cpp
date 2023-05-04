#include "./DictatNeuronPop.hpp"
#include "DictatNeuronPop.hpp"

DictatNeuronPop::DictatNeuronPop(GlobalSimInfo *info, int id) : NeuronPop (info, id)
{
    //Constructor
    inputFileAddress = info->pathTo_inputFile + "DictatNeuronPop_"+std::to_string(id)+"_spikers.txt";
    inputStreamFile.open(inputFileAddress, std::ifstream::in);
    //Then just open and read, with bools for raw spikers or instruction-based
}

void DictatNeuronPop::LoadParameters(std::vector<std::string> *input)
{
    NeuronPop::LoadParameters(input);//For now this is commented out to avoid confusion, and prevent from recording output in DictatNeuronPop

    std::string              name,token;
    std::vector<std::string> values;

    for(std::vector<std::string>::iterator it = (*input).begin(); it != (*input).end(); ++it) {

        SplitString(&(*it),&name,&values);

        if(name.find("inputFile") != std::string::npos){
            if(values.at(0).find("spiker") != std::string::npos){
                spikerListFiringBasedBool=true;
            } else if (values.at(0).find("instruction") != std::string::npos){
                instructionBasedBool=true;
                inputInstructions.resize(noNeurons);
                activeInstructions.resize(noNeurons);
                ReadInstructionsFromFile();
            } else {
                throw; //A non-existant type should throw
            }
        }else if(name.find("poissonLike") != std::string::npos){
            if(values.at(0).find("true") != std::string::npos){
                poissonLikeFiring=true;
                std::uniform_int_distribution<int> distribution(0,INT32_MAX);
                int seed = distribution(info->globalGenerator);
                randomGenerator = std::mt19937(seed);
            } else if (values.at(0).find("false") != std::string::npos){
                poissonLikeFiring=false;
            } else {
                throw; //A non-existant type should throw
            }
        }  
    }
}

void DictatNeuronPop::SaveParameters(std::ofstream *stream)
{

    std::string id = "neurons_" + std::to_string(GetId());

    NeuronPop::SaveParameters(stream);

    *stream <<  id + "_inputFile\t\t\t\t\t";
    if (spikerListFiringBasedBool){
        *stream<<"spiker";
    } else if (instructionBasedBool){
        *stream<<"instruction";
    } else {
        *stream<<"None";
    }
    *stream <<  "\t#Write 'instruction' for instruction-based, and 'spiker' for spiker-based \n";

    *stream <<  id + "_poissonLike\t\t\t\t\t";
    if (poissonLikeFiring){
        *stream<<"true";
    } else {
        *stream<<"false";
    }
    *stream <<  "\t#Write true vs false to indicate if the firing reproduced from instructions is poisson-like or periodic. \n";
}

void DictatNeuronPop::advect(std::vector<double> *synaptic_dV)
{
    ClearSpiker();
    if (spikerListFiringBasedBool){
        ReadSpikersFromFile();
    } else if(instructionBasedBool){
        if (poissonLikeFiring){
            GeneratePoissonSpikersFromInstructions();
        } else {
            GenerateRegularSpikersFromInstructions();
        }
    } else {
        throw;
    }
}

void DictatNeuronPop::ReadInstructionsFromFile()
{
    /*Instruction files for a population should be written in the following format:
    > neuronid starttime1 endtime1 fire_every_n_steps
    And ALWAYS time ordered with correct intervals (wrong starttimes will essentially change the first AP firing timestep or total phase)->(timeste-starttime)%frequency
    */
   //Runs only once in LP
   char*    entry = new char[2048];
    while(inputStreamFile.getline(entry, 256)){
        if(entry[0] == '#'){
            delete[] entry;
            entry = new char[2048];
            continue;
        }
        else if (entry[0] == '>') {
            FileEntry s_entry{ std::move(stringToFileEntry(std::move(static_cast<std::string>(entry)))) };
            inputInstructions.at(std::stoi(s_entry.values.at(0))).emplace_back(Instruction(std::stoi(s_entry.values.at(0)), std::stod(s_entry.values.at(1)), std::stod(s_entry.values.at(2)), std::stod(s_entry.values.at(3)), info->dt));
            delete[] entry;
            entry = new char[2048];
        } else if (entry[0] == ' ' || entry[0] == '\n') {
            delete[] entry;
            entry = new char[2048];
            continue;
        } else {
            std::cout<<"Reading error. The unexpected read input was: "<< entry <<"\n";
            throw;
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
    //         if ((!instruction.completed) && (instruction.startTimeStep>info->time_step)){ //Here not >= to avoid a spike of frequency when changing instructions
    //             if (((info->time_step-instruction.startTimeStep)%instruction.fireEveryNSteps)==0){
    //                 spiker.push_back(neuronId);
    //                 if (instruction.endTimeStep>=info->time_step){
    //                 instruction.completed=true;
    //             }
    //             break;
    //             } 
    //         }
    //     }
    // }
    //loop option 2 (efficient but probably worse in exception handling)
    for (int neuronId = 0; neuronId<static_cast<int>(activeInstructions.size()); neuronId++){
        Instruction& instruction = inputInstructions.at(neuronId).at(activeInstructions.at(neuronId));
        if (instruction.endTimeStep<=info->time_step){
            if (!(instruction.last)){
                activeInstructions.at(neuronId)++;
            }
            instruction.completed=true;
        }
        if (!instruction.off && ((info->time_step-instruction.startTimeStep)%instruction.fireEveryNSteps)==0){
            //Instruction.off condition is there to avoid doing a modulus with zero
            if (instruction.completed){continue;}
            if (instruction.startTimeStep<info->time_step) {spiker.push_back(neuronId);}
        }
    }
}

void DictatNeuronPop::GeneratePoissonSpikersFromInstructions()
{
        for (int neuronId = 0; neuronId<static_cast<int>(activeInstructions.size()); neuronId++){
        Instruction& instruction = inputInstructions.at(neuronId).at(activeInstructions.at(neuronId));
        if (instruction.endTimeStep<=info->time_step){
            if (!(instruction.last)){
                activeInstructions.at(neuronId)++;
            }
            instruction.completed=true;
        }
        if (uniformDistribution(randomGenerator)<instruction.lambdaPoisson){ //Here because there is no modulus, there is no need for checking the instruction.off condition
            if (instruction.completed){continue;}
            if (instruction.startTimeStep<info->time_step) {spiker.push_back(neuronId);}
        }
    }
}

void DictatNeuronPop::ReadSpikersFromFile()
{
    /*Spikers files for a population should be written in the following format:
    > time neuron1, neuron2, neuron 4.... (neurons that spike in time_step)
    *NEVER PUT COMMENTS OR EXTRA_LINES IN THIS TYPE OF FILE!
    *And be careful with your dt, as this will alter behaviour
    *This file should be generated automatically and tested against Test11
    *This function runs once per timestep
    */
   char* entry = new char[2048];
    inputStreamFile.getline(entry, 256);
        if(entry[0] == '>'){
            FileEntry s_entry {stringToFileEntry(std::move(static_cast<std::string>(entry)))};
            long timeStep{std::lround(std::stod(s_entry.values.at(0))/info->dt)};
            if (timeStep != info->time_step){
                std::cout<<"There was a reading alignment error";
                throw;
            } else {
                for (std::string value : s_entry.values){
                    spiker.push_back(std::stoi(value));
                }
            }
        } else {
            std::cout<<"Reading error. The unexpected read input was: "<< entry <<"\n";
            throw;
            }
}

Instruction::Instruction(int neuronId, double startTime, double endTime, double frequency, double dt): neuronId{neuronId}, startTimeStep{std::lround(startTime/dt)}, endTimeStep{std::lround(endTime/dt)}, frequency{frequency}, lambdaPoisson{frequency*dt}
{
    //Constructor
    if (frequency < std::numeric_limits<double>::epsilon()){ //Zero comparison to avoid division by zero
        this->off=true;
    } else {
        fireEveryNSteps=std::lround((1/frequency)/dt); //Conversion from frequency to timestep period.
        if (fireEveryNSteps==0){ //If the frequency is close
            std::cout<<"\n"<<"FATAL ERROR: YOU CHOSE A FREQUENCY THAT IS TOO HIGH FOR NEURON "<<std::to_string(neuronId)<<"\n\n\n"<<"**********************************";
            throw;
        }
    }
}
