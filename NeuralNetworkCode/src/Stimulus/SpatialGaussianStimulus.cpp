
//  Created by Pierre Ekelmans on 20/03/2018.
//  Copyright © 2018 Pierre Ekelmans. All rights reserved.
//

#include "SpatialGaussianStimulus.hpp"


SpatialGaussianStimulus::SpatialGaussianStimulus(NeuronPopSample *neur,std::vector<std::string> *input,GlobalSimInfo  * info):Stimulus(neur,info){

    //******************************
    //***** Default values *********
    //******************************
    seed         = 0;
    //******************************

    LoadParameters(input);
    generator  = std::mt19937(seed);
}


void SpatialGaussianStimulus::LoadParameters(std::vector<std::string> *input){

    //Stimulus::LoadParameters(input);
    std::string              name,token;
    std::vector<std::string> values;
    int                      P = neurons->GetTotalPopulations();
	int Gaussian_i;
    signal_step s;

	Ngauss = 1;
	maxCurrent.resize(Ngauss);
	sigmaCurrent_t.resize(Ngauss);
	sigmaCurrent_x.resize(Ngauss);
	GPos_X.resize(Ngauss);
	GPos_Y.resize(Ngauss);
	GPos_X[0] = 0.5;
	GPos_Y[0] = 0.5;

    for(std::vector<std::string>::iterator it = (*input).begin(); it != (*input).end(); ++it) {

        SplitString(&(*it),&name,&values);
		// By default, unless n is defined (see following input load), there is only one Gaussian, centered in 0.5 0.5
		if ((name.find("NumberOfGaussians") != std::string::npos)) {
			Ngauss = std::stoi(values.at(0));
/*			if (Ngauss > 9) {
				std::cout << "Only up to 9 Gaussians are supported at this point";
				Ngauss = 9;
			}*/
			maxCurrent.resize(Ngauss);
			sigmaCurrent_t.resize(Ngauss);
			sigmaCurrent_x.resize(Ngauss);
			GPos_X.resize(Ngauss);
			GPos_Y.resize(Ngauss);
		}

		else if ((name.find("X_position") != std::string::npos)) {
			for (int input = 0;input < Ngauss;input++) {
				GPos_X[input] = std::stod(values.at(input));
			}
		}

		else if ((name.find("Y_position") != std::string::npos) && (info->Dimensions==2)) {
			for (int input = 0;input < Ngauss;input++) {
				GPos_Y[input] = std::stod(values.at(input));
			}
		}

        else if((name.find("maxCurrent") != std::string::npos)){
			Gaussian_i = 1;
			if (name.find("t_") != std::string::npos) {
				Gaussian_i = std::stoi(name.substr(name.find("t_") + 2, std::string::npos));
			}
			s.values.resize(P);
			for (int input = 0;input < P;input++) {
				s.values[input] = std::stod(values.at(input));
			}
            if(is_double(values.at(P)))
				s.end_time = static_cast<int>(std::round(std::stod(values.at(P)) / info->dt));
            else
                s.end_time = INT_MAX;
			if ((s.end_time < 0) || (s.end_time >std::round(info->simulationTime / info->dt)))
				s.end_time = static_cast<int>(std::round(info->simulationTime / info->dt));
            maxCurrent[Gaussian_i-1].push_back(s);
        }

		else if ((name.find("sigmaCurrent_x") != std::string::npos)) {
			Gaussian_i = 1;
			if (name.find("x_") != std::string::npos) {
				Gaussian_i = std::stoi(name.substr(name.find("x_") + 2, std::string::npos));
			}
			s.values.resize(1);
			s.values[0] = std::stod(values.at(0));
			if (is_double(values.at(1)))
				s.end_time = static_cast<int>(std::round(std::stod(values.at(1)) / info->dt));
			else
				s.end_time = INT_MAX;
			if ((s.end_time < 0) || (s.end_time >std::round(info->simulationTime / info->dt)))
				s.end_time = static_cast<int>(std::round(info->simulationTime / info->dt));
			sigmaCurrent_x[Gaussian_i-1].push_back(s);
		}

        else if((name.find("sigmaCurrent_t") != std::string::npos)){
			Gaussian_i = 1;
			if (name.find("_t_") != std::string::npos) {
				Gaussian_i = std::stoi(name.substr(name.find("_t_") + 3, std::string::npos));
			}
			s.values.resize(P);
			for (int input = 0;input < P;input++) {
				s.values[input] = std::stod(values.at(input));
			}
			if (is_double(values.at(P)))
				s.end_time = static_cast<int>(std::round(std::stod(values.at(P)) / info->dt));
			else
				s.end_time = INT_MAX;
			if ((s.end_time < 0) || (s.end_time >std::round(info->simulationTime / info->dt)))
				s.end_time = static_cast<int>(std::round(info->simulationTime / info->dt));
			sigmaCurrent_t[Gaussian_i-1].push_back(s);
        }
		else if ((name.find("Background_Noise") != std::string::npos)) {
			s.values.resize(P);
			for (int input = 0;input < P;input++)
				s.values[input] = std::stod(values.at(input));
			if (is_double(values.at(P)))
				s.end_time = static_cast<int>(std::round(std::stod(values.at(P)) / info->dt));
			else
				s.end_time = INT_MAX;
			if ((s.end_time < 0) || (s.end_time >std::round(info->simulationTime / info->dt)))
				s.end_time = static_cast<int>(std::round(info->simulationTime / info->dt));
			backgroundNoise.push_back(s);
		}
    }

    if(info->globalSeed != -1){
        std::uniform_int_distribution<int> distribution(0,INT32_MAX);
        seed = distribution(info->globalGenerator);
    }
	else {
		std::cout << "Please use a positive global seed to use SpatialGaussianStimulus";
		throw "Someone attempted to use SpatialGaussianStimulus without a global seed";
	}
}


void SpatialGaussianStimulus::SaveParameters(std::ofstream * stream){

    int P        = neurons->GetTotalPopulations();
    Stimulus::SaveParameters(stream);

	*stream << "stimulus_NumberOfGaussians           ";
	*stream << std::to_string(Ngauss) << "\n";

	*stream << "stimulus_X_position                  ";
	for (int Gauss_index = 0; Gauss_index < Ngauss;Gauss_index++)
		*stream << std::to_string(GPos_X[Gauss_index]) << "\t ";
	*stream << "Position of each Gaussian on the X axis (between 0 and 1)\n";

	if (info->Dimensions == 2) {
		*stream << "stimulus_Y_position                  ";
		for (int Gauss_index = 0; Gauss_index < Ngauss;Gauss_index++)
			*stream << std::to_string(GPos_Y[Gauss_index]) << "\t ";
		*stream << "Position of each Gaussian on the Y axis (between 0 and 1)\n";
	}

	for (int Gauss_index = 0; Gauss_index < Ngauss;Gauss_index++) {

		for (auto const &s : maxCurrent[Gauss_index]) {
			*stream << "stimulus_maxCurrent_" << std::to_string(Gauss_index + 1) << "                ";
			for (int i = 0;i < P;i++)
				*stream << std::to_string(s.values[i]) << "\t ";
			*stream << std::to_string(static_cast<double>(s.end_time)*info->dt) << " \t";
			*stream << " [column i: input to neurons of population i, at the center of the gaussian, last column: time until which input is set. Dimensions: [mV/sec , secs.]\n";
		}

		for (auto const &s : sigmaCurrent_t[Gauss_index]) {
			*stream << "stimulus_sigmaCurrent_t_" << std::to_string(Gauss_index + 1) << "            ";
			for (int i = 0;i < P;i++)
				*stream << std::to_string(s.values[i]) << "\t ";
			*stream << std::to_string(static_cast<double>(s.end_time)*info->dt) << " \t";
			*stream << " [column i: relative input noise to population i (relative to the mean current), last column: time until which input is set. Dimensions: [ -  , secs.]\n";
		}

		for (auto const &s : sigmaCurrent_x[Gauss_index]) {
			*stream << "stimulus_sigmaCurrent_x_" << std::to_string(Gauss_index + 1) << "            ";
			*stream << std::to_string(s.values[0]) << "\t ";//the width of the input is the same for all popluations
			*stream << std::to_string(static_cast<double>(s.end_time)*info->dt) << " \t";
			*stream << " [column 1: spatial spread (std of the Gaussian) of the input to all populations, last column: time until which input is set. Dimensions: [mm , secs.]\n";
		}
	}

	for (auto const &s : backgroundNoise) {
		*stream << "stimulus_Background_Noise            ";
		for (int i = 0;i<P;i++)
		*stream << std::to_string(s.values[i]) << "\t ";
	*stream << std::to_string(static_cast<double>(s.end_time)*info->dt) << " \t";
	*stream << "Noise applied in the whole domain [mV/sqrt(sec) , secs.]\n";
	}

	*stream <<  "#\t\tRI_{i,ext}/tau_m*dt = meanCurrent_i*dt*exp(-d{i}^2/(2sigmaCurrent_x)) + sqrt(dt)*sigmaCurrent_t_i*NormalDistribution(0,1)\t Where d{i} is the distance of neuron i to the center of the domain\n";
}


void SpatialGaussianStimulus::SetSignalArray(){

	double mu_max, sigmax, sigmat;
	double s;
    double dt      = info->dt;
    int    P       = neurons->GetTotalPopulations();
    double sqrt_dt = sqrt(dt);
    long   t_step  = info->time_step;
	double   mu;

	double d;
    std::normal_distribution<double> distribution(0,1);
	signal_step   *max_step_current;
	signal_step   *sigma_step_current_x;
	signal_step   *sigma_step_current_t;
	signal_step	*background_step;

	background_step = &backgroundNoise.at(0);
	while ((t_step > background_step->end_time) && (background_step != &backgroundNoise.back()))
		background_step++;

	for (int Gauss_index = 0; Gauss_index < Ngauss;Gauss_index++) {
		max_step_current = &maxCurrent[Gauss_index].at(0);
		sigma_step_current_x = &sigmaCurrent_x[Gauss_index].at(0);
		sigma_step_current_t = &sigmaCurrent_t[Gauss_index].at(0);

		while ((t_step > max_step_current->end_time) && (max_step_current != &maxCurrent[Gauss_index].back()))
			max_step_current++;

		while ((t_step > sigma_step_current_t->end_time) && (sigma_step_current_t != &sigmaCurrent_t[Gauss_index].back()))
			sigma_step_current_t++;

		while ((t_step > sigma_step_current_x->end_time) && (sigma_step_current_x != &sigmaCurrent_x[Gauss_index].back()))
			sigma_step_current_x++;

		sigmax = sigma_step_current_x->values.at(0);
		for (int pop = 0;pop<P;pop++) {
			mu_max = max_step_current->values.at(pop);
			sigmat = sigma_step_current_t->values.at(pop);
			s = GetScaling(pop);

			for (unsigned int i = 0;i < neurons->GetNeuronsPop(pop); i++) {
				d = pow(neurons->GetX_Pos(pop, i) - GPos_X[Gauss_index]*info->Lx, 2) + pow(neurons->GetY_Pos(pop, i) - GPos_Y[Gauss_index]*info->Ly, 2);
				mu = mu_max * dt*pow(s, -(info->networkScaling_synStrength))*exp(-d / (2 * pow(sigmax, 2)));
				//mu = 0;//For repeating boundaries
				//for (int x_rep = -1;x_rep <= 1;x_rep++) {
				//	for (int y_rep = -1;y_rep <= 1;y_rep++) {
				//		d = pow(neurons->GetX_Pos(pop, i) - info->Lx *(GPos_X[Gauss_index] + x_rep), 2) + pow(neurons->GetY_Pos(pop, i) - info->Ly *(GPos_Y[Gauss_index] + y_rep), 2);
				//		if (info->Dimensions == 2 || y_rep == 0)
				//			mu = mu + mu_max * dt*pow(s, -(info->networkScaling_synStrength))*exp(-d / (2 * pow(sigmax, 2)));
				//	}
				//}
				if (Gauss_index == 0)
					signal_array[pop][i] = mu * (1 + sigmat * distribution(generator)) + background_step->values.at(pop) * sqrt_dt * distribution(generator);
				else
					signal_array[pop][i] += mu* (1 + sigmat * distribution(generator));
			}
		}

	}
}

void SpatialGaussianStimulus::Update(std::vector<std::vector<double>> * synaptic_dV)
{
    SetSignalArray();
    Stimulus::Update(synaptic_dV);
}
