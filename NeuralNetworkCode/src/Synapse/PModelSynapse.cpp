#include "PModelSynapse.hpp"

PModelSynapse::PModelSynapse(PopPtr targetPop, PopPtr sourcePop, GlobalSimInfo *infoGlobal) : Synapse(targetPop, sourcePop, infoGlobal) {
  hasPlasticityModel = true;
}
std::vector<double> PModelSynapse::AdvectSpikers(NeuronInt spiker) {
  // This is basically the CurrentSynapse AdvectSpikers, but for synapse spine pointers not nullptr
  NeuronInt           noTargets{GetNoTargetedSynapses(spiker)};
  std::vector<double> currents(noTargets, 0.0);
  for (NeuronInt targetNeuronIndex : std::ranges::views::iota(0, noTargets)) {
    // for(NeuronInt targetNeuronIndex{0};targetNeuronIndex<noTargets;targetNeuronIndex++){
    this->targetPop->RecordExcitatorySynapticSpike(targetSpineList.at(spiker).at(targetNeuronIndex).first,
                                                   targetSpineList.at(spiker).at(targetNeuronIndex).second);
    currents.at(targetNeuronIndex) = targetSpineList.at(spiker).at(targetNeuronIndex).second->GetWeight();
  }
  return currents;
}

void PModelSynapse::SaveParameters(std::ofstream &wParameterStream, std::string idString) const {
  Synapse::SaveParameters(wParameterStream, idString);

  if (this->targetPop->HasPlasticityModel() && !branchTarget.setOfPositions.empty()) {
    wParameterStream << idString << "positionsSet\t\t\t\t\t";
    for (int position : branchTarget.setOfPositions) {
      wParameterStream << std::to_string(position) << ' ';
    }
    wParameterStream << "\t"
                     << "# Positions to place synapses";
  }
  // Synaptic targeting
  if (this->targetPop->IsBranchedBool()) {
    wParameterStream << idString << "targetBranch\t\t\t\t\t";
    if (this->branchTarget.randomTargetBranch) {
      wParameterStream << "random"; // Missing comments on what this is supposed to do
    } else if (this->branchTarget.setTargetBranch) {
      wParameterStream << std::to_string(this->branchTarget.targetBranch); // Missing comments on what
                                                                           // this is supposed to do
    } else if (this->branchTarget.orderedTargetBranch) {
      wParameterStream << "ordered";
    } else {
      wParameterStream << "none"; // Missing comments on what this is supposed to do
    }
    wParameterStream << "\t\t\t#You can target branches in an 'ordered' "
                        "manner (0,1,2...), "
                        "'random', or set (if you input a number). Put "
                        "none if the HS does "
                        "not used branched morphology\n";
    wParameterStream << idString + "pmodel_"
                     << "slotOrder\t\t";
    if (this->branchTarget.firstSlotTrueLastSlotFalse) {
      wParameterStream << "first\t";
    } else {
      wParameterStream << "last\t";
    }
    wParameterStream << "\t"
                     << "#'first' synapse allocation will allocate synapses from the "
                        "beggining to the end of the available slots. 'last' will do "
                        "the "
                        "opposite. This only makes sense in ordered allocation\n";
  }
  targetPop->SavePlasticityModel(wParameterStream, idString + "pmodel_"); //!!!!

  targetPop->ModelSaved();
  // if (geometry != nullptr){
  //     wParameterStream << "########### Connectivity ################\n";
  //     geometry->SaveParameters(wParameterStream,idString);
  // }
}

void PModelSynapse::LoadParameters(const std::vector<FileEntry> &hybridParameters) {
  branchTarget.setOfPositions.clear();
  for (auto &[parameterName, parameterValues] : hybridParameters) {
    if (parameterName.find("targetBranch") != std::string::npos) {
      if (parameterValues.at(0).find("random") != std::string::npos) {
        this->branchTarget.randomTargetBranch = true;
      } else if (parameterValues.at(0).find("none") != std::string::npos) {
        // Here nothing is done to handle the case where we do not used
        // branched Checking if the input is an integer would take more
        // time than checking if None is in a single string with 5
        // characters max (I think).
      } else if (parameterValues.at(0).find("ordered") != std::string::npos) {
        this->branchTarget.orderedTargetBranch = true;
      } else {
        this->branchTarget.targetBranch = std::stoi(parameterValues.at(0));
        if (this->branchTarget.targetBranch >= targetPop->GetNoBranches()) {
          throw "Targeted branch does not exist";
        }
        this->branchTarget.setTargetBranch = true;
        // Missing exception management for when the input is not an
        // integer.
      }
    } else if (parameterName.find("subregion") != std::string::npos) {
      this->branchTarget.DendriticSubRegion = parameterValues.at(0).at(0); // Re-think this char array wacky stuff
    } else if (parameterName.find("positionsSet") != std::string::npos && !parameterValues.empty()) {
      for (std::string param : parameterValues) {
        branchTarget.setOfPositions.push_back(std::stoi(param));
      }
    } else if (parameterName.find("slotOrder") != std::string::npos) {
      if (parameterValues.at(0).find("first") != std::string::npos) {
        this->branchTarget.firstSlotTrueLastSlotFalse = true;
      } else if (parameterValues.at(0).find("last") != std::string::npos) {
        this->branchTarget.firstSlotTrueLastSlotFalse = false;
      }
    }
  }
  targetPop->LoadPlasticityModel(FilterStringEntries(hybridParameters, "pmodel"));
  this->ignoreJDParameters = targetPop->ignoreJDParameters();
  Synapse::LoadParameters(hybridParameters);
  Dmax = targetPop->GetMaxGapDelay(delayPerMicrometer);
}
void PModelSynapse::AllocateSynapse(NeuronInt targetNeuron, NeuronInt sourceNeuron) {
  AllocateSynapseWithPlasticity(targetNeuron, sourceNeuron);
}
std::string PModelSynapse::GetDataHeader(int dataColumn) {
  return "#" + std::to_string(dataColumn) + " J_" + GetIdStr() + " (mV)\n";
}

std::string PModelSynapse::GetUnhashedDataHeader() const {
  return "J_" + GetIdStr() + "\t";
}

std::vector<double> PModelSynapse::GetSynapticState(NeuronInt sourceNeuron) const {
  std::vector<double> value(1);
  double              Jsum = 0;
  // get average coupling strength
  Jsum        = std::accumulate(targetSpineList.at(sourceNeuron).begin(), targetSpineList.at(sourceNeuron).end(), 0.0,
                                [](double Jsum, std::pair<NeuronInt, BaseSpinePtr> synapse) { return Jsum + synapse.second->GetWeight(); });
  value.at(0) = Jsum / static_cast<double>(this->GetNoTargetedSynapses(sourceNeuron));
  // value[0] = GetCouplingStrength()*static_cast<double>(this->GetNumberOfPostsynapticTargets(pre_neuron));
  return value;
}
