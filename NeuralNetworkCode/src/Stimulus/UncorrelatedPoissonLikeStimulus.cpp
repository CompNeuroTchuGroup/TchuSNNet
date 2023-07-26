#include "UncorrelatedPoissonLikeStimulus.hpp"

/* When the stimulus has to be changed, the stimulus value at the end
 * is deleted such that at the end of nextStimStep there is the current
 * stimulus value.
 * The last value of nextStimTimeStep is also deleted such that
 * at the end there is the timestep for the next stimulus change.
 */

UncorrelatedPoissonLikeStimulus::UncorrelatedPoissonLikeStimulus(std::shared_ptr<NeuronPopSample> neurons,std::vector<FileEntry>& stimulusParameters, GlobalSimInfo*  infoGlobal): Stimulus(neurons,infoGlobal){

    PopInt totalNeuronPops             = neurons->GetTotalPopulations();
    // tableEntries     = 10;
    J_External.resize(totalNeuronPops, 0.0);
    externalCurrents.resize(totalNeuronPops, 0.0);
    // seed              = -1;

    // for(PopInt neuronPop : std::ranges::views::iota(0, totalNeuronPops))
    //     J_External[i] = 0.0;

    LoadParameters(stimulusParameters);
}

void UncorrelatedPoissonLikeStimulus::LoadParameters(const std::vector<FileEntry>& stimulusParameters){

    Stimulus::LoadParameters(stimulusParameters);
    PopInt  totalNeuronPops = neurons->GetTotalPopulations();

    for(auto& [parameterName, parameterValues] : stimulusParameters) {

        if((parameterName.find("virtualExternalNeurons") != std::string::npos) ||
           (parameterName.find("noExternalNeurons") != std::string::npos)){
            this->noExternalNeurons = std::stoi(parameterValues.at(0));
        // } else if(parameterName.find("seed") != std::string::npos){
        //     SetSeed(static_cast<int>(std::stod(parameterValues.at(0))));
        } else if(parameterName.find("stimulus_step") != std::string::npos){
            AddStimulusStep(std::stod(parameterValues.at(0)), std::stod(parameterValues.at(1)));
        // } else if((parameterName.find("Poisson_table_entries") != std::string::npos) ||
        //         (parameterName.find("PoissonTableEntries") != std::string::npos)){
        //     this->tableEntries = static_cast<int>(std::stod(parameterValues.at(0)));
		} else if (parameterName.find("J_X") != std::string::npos) {
			for (PopInt neuronPop : std::ranges::views::iota(0,totalNeuronPops)){
				try{
					this->J_External.at(neuronPop) = std::stod(parameterValues.at(neuronPop));
				} catch (...){
					throw "J_X parameters do not account for all populations";
				}
			}
        }
    }
	if(infoGlobal->isMock){
		return;
	}
    PostLoadParameters();
}
void UncorrelatedPoissonLikeStimulus::SaveParameters(std::ofstream& wParameterStream) const{

    PopInt totalNeuronPops        = neurons->GetTotalPopulations();
    Stimulus::SaveParameters(wParameterStream);

    wParameterStream <<  "stimulus_noExternalNeurons           " << std::to_string(noExternalNeurons)  << "\n";
    // if(infoGlobal->globalSeed == -1){
    //     wParameterStream <<  "stimulus_seed                        " << std::to_string(seed)  << "\n";
    // }
    // wParameterStream <<  "stimulus_PoissonTableEntries         " << std::to_string(GetTableEntries())  << "\n";
    wParameterStream <<  "stimulus_J_X                         ";
    for (PopInt neuronPop : std::ranges::views::iota(0,totalNeuronPops)){
        wParameterStream << std::to_string(J_External.at(neuronPop)) << "\t";
    }
    wParameterStream << " #dmV/Spike\n";

    for(int step : std::ranges::views::iota(0,GetStimulusNoSteps())){
        wParameterStream <<  "stimulus_step                        " << std::to_string(GetStimulusStepEndTime(step)*infoGlobal->dtTimestep)  << "\t" << std::to_string(GetStimulusStep(step))  << "\t [t (secs.) -- Hz]\n";
    }
    wParameterStream <<  "#\t\t" << IDstringUncorrelatedStimulus << ": noExternalNeurons external neurons with (poisson) firing rate defined in stimulus_step are connected with population X with strength J_X.\n";
}

void UncorrelatedPoissonLikeStimulus::PostLoadParameters(){
	std::cout<<"\nSetting up the stimulus class...";
    if(stimulusSteps.empty()){ //Nothing should happen if there is no structs in the vector
    	wellDefined=false;
        std::cout<<"Stimulus class was ill-defined\n"<<std::endl;
        return;
    }
    std::sort(stimulusSteps.begin(), stimulusSteps.end(), [](StepStruct& step1,StepStruct& step2){
		return step1.endTimeStep<step2.endTimeStep;
	});
    currentStep=stimulusSteps.data();
	for(PopInt neuronPop : std::ranges::views::iota(0,neurons->GetTotalPopulations())){
		externalCurrents.at(neuronPop)= J_External.at(neuronPop) * pow(GetScaling(neuronPop),(infoGlobal->networkScaling_synStrength));
	}
    //Previously in SetTableEntries
    //double signal = infoGlobal->dtTimestep*static_cast<double>(noExternalNeurons)*stimulusSteps.front().parameterValues.at(0);

    poissonDistr=std::poisson_distribution<int>(infoGlobal->dtTimestep*static_cast<double>(noExternalNeurons)*currentStep->parameterValues.at(0));
    std::cout<<"Done.\n"<<std::endl;
}

void UncorrelatedPoissonLikeStimulus::AddStimulusStep(double endTime,double stimStep){
    StepStruct step;
	step.endTimeStep=static_cast<TStepInt>(std::lround(endTime/infoGlobal->dtTimestep));
	step.parameterValues.push_back(stimStep);
	stimulusSteps.push_back(step);
    // TStepInt	endTimeStep{ static_cast<TStepInt>(std::lround(endTime / infoGlobal->dtTimestep)) };
    // if(nextStimStep.empty())
    // {
    //   nextStimTimeStep.push_back(endTimeStep);
    //   nextStimStep.push_back(stimStep);
    // }
    // else
    // {
    //   if(endTimeStep > nextStimTimeStep.back())
    //   {
    //     std::vector<double> temp_time_step(nextStimTimeStep);
    //     std::vector<double> temp_stimulus(nextStimStep);

    //     nextStimTimeStep.clear();
    //     nextStimStep.clear();

    //     bool done_flag = false;
    //     for(signed i = 0; i < temp_time_step.size(); i++)
    //     {
    //       if(!done_flag && (endTimeStep > temp_time_step[i]))
    //       {
    //         nextStimTimeStep.push_back(endTimeStep);
    //         nextStimStep.push_back(stimStep);
    //         done_flag = true;
    //       }
    //       nextStimTimeStep.push_back(temp_time_step[i]);
    //       nextStimStep.push_back(temp_stimulus[i]);
    //     }
    //     if(nextStimTimeStep.size() == temp_time_step.size())
    //     {
    //       nextStimTimeStep.push_back(endTimeStep);
    //       nextStimStep.push_back(stimStep);
    //     }
    //   }
    //   else
    //   {
    //     nextStimTimeStep.push_back(endTimeStep);
    //     nextStimStep.push_back(stimStep);
    //   }
    // }
}

void UncorrelatedPoissonLikeStimulus::UpdatePoissonTable() {
   	// StepStruct& step = *std::find_if(stimulusSteps.begin(), stimulusSteps.end(), [this](StepStruct step){
	// 	//Confusing unless you think about a while loop with the condition
	// 	return (infoGlobal->timeStep<=step.endTimeStep) || (step.lastStep);
	// });
    if((currentStep != &stimulusSteps.back()) && (currentStep->endTimeStep<=infoGlobal->dtTimestep)){
		currentStep++;
    	poissonDistr=std::poisson_distribution<int>(infoGlobal->dtTimestep*static_cast<double>(noExternalNeurons)*currentStep->parameterValues.at(0));
	}
    // double  signal;
    // double  dtTimestep = infoGlobal->dtTimestep;
    // if(!nextStimTimeStep.empty())
    // {
    //    while(nextStimTimeStep.back() <= infoGlobal->timeStep)
    //     {
    //         nextStimStep.pop_back();
    //         nextStimTimeStep.pop_back();
    //         signal = dtTimestep*static_cast<double>(noExternalNeurons)*nextStimStep.back();
    //         FillPoissonValueTable(signal);

    //         if(nextStimTimeStep.empty())
    //             break;
    //     }
    // }
}
void UncorrelatedPoissonLikeStimulus::SetSignalMatrix() {

    UpdatePoissonTable();

    for(PopInt neuronPop : std::ranges::views::iota(0,neurons->GetTotalPopulations())){
        for(NeuronInt neuron : std::ranges::views::iota(0,neurons->GetNeuronsPop(neuronPop))){
            signalMatrix.at(neuronPop).at(neuron)   = externalCurrents.at(neuronPop)*poissonDistr(generator);
        }
    }
}

void UncorrelatedPoissonLikeStimulus::Update(std::vector<std::vector<double>>& synaptic_dV){
    if(!wellDefined){
        return;
    }
    SetSignalMatrix();
    Stimulus::Update(synaptic_dV);
}

double UncorrelatedPoissonLikeStimulus::GetScaling(PopInt neuronPop) const {
        double avgConnExtPop;//Average connections of the external populations
        if (infoGlobal->networkScaling_mode == 0){
            avgConnExtPop = static_cast<double>(noExternalNeurons);
        } else if (infoGlobal->networkScaling_mode == 1){
            avgConnExtPop = static_cast<double>(neurons->GetTotalNeurons());
        } else{
            throw "ERROR: GetExternalCouplingStrength";
        }
        return avgConnExtPop;
}
// fills the poissonValueTable that is a tabular version of the poisson distribution
// with the needed mean that is firing_rate*dt*number_of_neurons
// inline void UncorrelatedPoissonLikeStimulus::FillPoissonValueTable(double mu)
// {

// 	int value                    = 0;
// 	double probability           = exp(-mu);
// 	double cumulated_probability = exp(-mu);

// 	for(int i = 0; i < tableEntries; ++i)
// 	{
// 		if(cumulated_probability < static_cast<double>(i)/static_cast<double>(tableEntries))
// 		{
// 			++value;
// 			probability            = probability*mu/static_cast<double>(value);
// 			cumulated_probability += probability;
// 		}
// 		poissonValueTable[i] = static_cast<double>(value);
// 	}
// }



