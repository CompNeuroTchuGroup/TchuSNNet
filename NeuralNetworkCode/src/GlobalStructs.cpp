#include "GlobalStructs.hpp"

DendriticSubRegion::DendriticSubRegion(char regionID, std::vector<int> branchesInRegion):
    regionID { regionID }, branchesInRegion { branchesInRegion } {
  // This is still under work
}

void FileEntry::RemoveCommentsInValues(char commentCharacter) {
  std::vector<std::string> uncommentedValues;
  for (std::string& storedValue : parameterValues) {
    if (storedValue.find(commentCharacter) != std::string::npos) {
      break;
    } else {
      uncommentedValues.push_back(storedValue);
    }
  }
  this->parameterValues = uncommentedValues;
}

Instruction::Instruction(FileEntry inputEntry, double dtTimestep):
    neuronId { std::stoi(inputEntry.parameterValues.at(0)) },
    startTimeStep { std::lround(std::stod(inputEntry.parameterValues.at(1)) / dtTimestep) },
    endTimeStep { std::lround(std::stod(inputEntry.parameterValues.at(2)) / dtTimestep) },
    frequency { std::stod(inputEntry.parameterValues.at(3)) }, firingProbability { frequency * dtTimestep } {
  // Constructor
  if (frequency < std::numeric_limits<double>::epsilon()) {  // Zero comparison to avoid division by zero
    this->off = true;
  } else {
    fireEveryNSteps = std::lround((1 / frequency) / dtTimestep);  // Conversion from frequency to timestep period.
    if (fireEveryNSteps == 0) {                                   // If the frequency is close
      std::cout << "\n"
                << "EXCEPTION: YOU CHOSE A FREQUENCY THAT IS TOO HIGH FOR NEURON " << std::to_string(neuronId) << "\n\n\n"
                << "**********************************";
      throw "EXCEPTION: YOU CHOSE A FREQUENCY THAT IS TOO HIGH FOR CURRENT DT IN DICTAT INPUT FILE";
    }
  }
}

// void ExponentialFit::fitExponential(std::vector<double>& x_data, std::vector<double>& y_data) {
//   int    n = x_data.size();
//   double delta_a, delta_b;
//   double error, error_new;

//   for (int iter = 0; iter < max_iterations; ++iter) {
//     double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0;

//     for (int i = 0; i < n; ++i) {
//       double y_fit    = exponential_func(x_data[i]);
//       double residual = y_data[i] - y_fit;

//       jacobian_exponential(x_data[i]);

//       sum1 += residual * dfda * dfda;
//       sum2 += residual * dfda * dfdb;
//       sum3 += residual * dfdb * dfda;
//       sum4 += residual * dfdb * dfdb;
//     }

//     delta_a = (sum4 * sum1 - sum2 * sum3) / (n * sum1 - sum2 * sum3);
//     delta_b = (n * sum2 - sum1 * sum3) / (n * sum1 - sum2 * sum3);

//     init_val += delta_a;
//     exp_rate += delta_b;

//     // Check convergence
//     error_new = 0.0;
//     for (int i = 0; i < n; ++i) {
//       double y_fit = exponential_func(x_data[i]);
//       error_new += (y_data[i] - y_fit) * (y_data[i] - y_fit);
//     }

//     if (std::fabs(error_new - error) < tolerance) {
//       break;
//     }

//     error = error_new;
//   }
// }

// void ExponentialFit::jacobian_exponential(double x) {
//   dfda = std::exp(exp_rate * x);
//   dfdb = init_val * x * std::exp(exp_rate * x);
// }
