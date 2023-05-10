
#include "WhiteNoiseLinear.hpp"


WhiteNoiseLinear::WhiteNoiseLinear(NeuronPopSample *neur, std::vector<std::string> *input, GlobalSimInfo  * info) :Stimulus(neur, info) {

	//******************************
	//***** Default values *********
	//******************************
	seed = 0;
	//******************************

	LoadParameters(input);
	generator = std::mt19937(seed);
}


void WhiteNoiseLinear::LoadParameters(std::vector<std::string> *input) {

	//Stimulus::LoadParameters(input);
	std::string              name, token;
	std::vector<std::string> values;
	unsigned int                      P = neurons->GetTotalPopulations();
	Step s;

	for (std::vector<std::string>::iterator it = (*input).begin(); it != (*input).end(); ++it) {

		SplitString(&(*it), &name, &values);

		if ((name.find("seed") != std::string::npos)) {
			seed = static_cast<unsigned int>(std::stod(values.at(0)));
		}
		else if ((name.find("meanCurrent") != std::string::npos)) {
			s.values.clear();
			for (unsigned int i = 0;i < P;i++) {
				s.values.push_back(std::stod(values.at(2 * static_cast<long long>(i))));
				s.values.push_back(std::stod(values.at(2 * static_cast<long long>(i) + 1)));
			}
			s.start_time = static_cast<int>(std::round(std::stod(values.at(2 * P)) / info->dt));
			s.end_time = static_cast<int>(std::round(std::stod(values.at(2 * P + 1)) / info->dt));

			meanCurrent.push_back(s);
		}
		else if ((name.find("sigmaCurrent") != std::string::npos)) {
			s.values.clear();
			for (unsigned int i = 0;i< P;i++){
				s.values.push_back(std::stod(values.at(2*i)));
				s.values.push_back(std::stod(values.at(2*i+1)));
				}
			s.start_time = static_cast<int>(std::round(std::stod(values.at(2 * P)) / info->dt));
			s.end_time = static_cast<int>(std::round(std::stod(values.at(2 * P + 1)) / info->dt));
			sigmaCurrent.push_back(s);
		}
	}

	if (info->globalSeed != -1) {
		std::uniform_int_distribution<int> distribution(0, INT32_MAX);
		seed = distribution(info->globalGenerator);
	}
}


void WhiteNoiseLinear::SaveParameters(std::ofstream * stream) {

	unsigned int P = neurons->GetTotalPopulations();
	Stimulus::SaveParameters(stream);

	if (info->globalSeed == -1) {
		*stream << "stimulus_seed                        " << std::to_string(seed) << "\n";
	}

	for (auto const &s : meanCurrent) {
		*stream << "stimulus_meanCurrent                 ";
		for (unsigned int i = 0;i < P;i++) {
			*stream << std::to_string(s.values.at(2 * i)) << "\t ";
			*stream << std::to_string(s.values.at(2 * i + 1)) << "\t ";
		}
		*stream << std::to_string(static_cast<double>(s.start_time)*info->dt) << " \t";
		*stream << std::to_string(static_cast<double>(s.end_time)*info->dt) << " \t";
		*stream << " [col 1: input to pop0 at t_0, col 2: pop1 at t_0, ... colP+1: pop1 t_f, ... col2P: popN t_f, t0, tf. Dimensions: [mV/sec , secs.]\n";
	}

	for (auto const &s : sigmaCurrent) {
		*stream << "stimulus_sigmaCurrent                ";
		for (unsigned int i = 0;i < P;i++) {
			*stream << std::to_string(s.values.at(2 * i)) << "\t ";
			*stream << std::to_string(s.values.at(2 * i + 1)) << "\t ";
		}
		*stream << std::to_string(static_cast<double>(s.start_time)*info->dt) << " \t";
		*stream << std::to_string(static_cast<double>(s.end_time)*info->dt) << " \t";
		*stream << " [col 1: input to pop0 at t_0, col 2: pop1 at t_0, ... colP+1: pop1 t_f, ... col2P: popN t_f, t0, tf. Dimensions: [mV/sqrt(sec) , secs.]\n";
	}

	*stream << "#\t\tRI_{i,ext}/tau_m*dt = meanCurrent_i*dt + sqrt(dt)*sigmaCurrent_i*NormalDistribution(0,1)\n";
}


void WhiteNoiseLinear::SetSignalArray() {

	double mean, sigma, s;
	double t0, tf;
	double mu0, muf, sig0,sigf;
	double dt = info->dt;
	unsigned int    P = neurons->GetTotalPopulations();
	double sqrt_dt = sqrt(dt);
	long   t_step = info->time_step;
	std::normal_distribution<double> distribution(0, 1);
	Step   *mean_step_current = &meanCurrent.at(0);
	Step   *sigma_step_current = &sigmaCurrent.at(0);;

	while ((t_step > mean_step_current->end_time) && (mean_step_current != &meanCurrent.back()))
		mean_step_current++;

	while ((t_step > sigma_step_current->end_time) && (sigma_step_current != &sigmaCurrent.back()))
		sigma_step_current++;

	for (unsigned int pop = 0;pop < P;pop++) {
		mu0 = mean_step_current->values.at(pop);
		muf = mean_step_current->values.at(pop + P);
		t0 = mean_step_current->start_time;
		tf = mean_step_current->end_time;
		mean = mu0 + (muf - mu0)*(t_step - t0) / (tf - t0);

		sig0 = sigma_step_current->values.at(pop);
		sigf = sigma_step_current->values.at(pop + P);
		t0 = sigma_step_current->start_time;
		tf = sigma_step_current->end_time;
		sigma = sig0 + (sigf - sig0)*(t_step - t0) / (tf - t0);
		s = GetScaling(pop);
		for (unsigned int i = 0;i < neurons->GetNeuronsPop(pop); i++){
			signal_array[pop][i] = mean * dt*pow(s, -(info->networkScaling_synStrength)) + sqrt_dt * sigma * distribution(generator);
		}
	}
}

void WhiteNoiseLinear::Update(std::vector<std::vector<double>> * synaptic_dV)
{
	SetSignalArray();
	Stimulus::Update(synaptic_dV);

}
