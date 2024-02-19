#include "./NoStimulus.hpp"
#include "NoStimulus.hpp"

void NoStimulus::SetSignalMatrix() {
}

void NoStimulus::PostLoadParameters() {
}

NoStimulus::NoStimulus(std::shared_ptr<NeuronPopSample> neurons, GlobalSimInfo *infoGlobal) : Stimulus(neurons, infoGlobal) {
}

void NoStimulus::Update(std::vector<std::vector<double>> &synaptic_dV) {
}

void NoStimulus::SaveParameters(std::ofstream &wParameterStream) const {
  Stimulus::SaveParameters(wParameterStream);
}

void NoStimulus::LoadParameters(const std::vector<FileEntry> &input) {
}
