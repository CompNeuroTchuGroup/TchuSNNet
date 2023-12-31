
#include "Synapse.hpp"
Synapse::Synapse(PopPtr targetPop, PopPtr sourcePop, GlobalSimInfo* infoGlobal) : sourcePop{sourcePop}, targetPop{targetPop},sourcePopID{sourcePop->GetId()}, targetPopID{targetPop->GetId()},infoGlobal{infoGlobal} {
    std::uniform_int_distribution<int> distribution(0,INT32_MAX);
    SetSeed(distribution(infoGlobal->globalGenerator));
    //********************************
    //***** Set Default parameterValues *******
    //********************************
    //********************************
    //********************************
    targetSpineList = std::vector<std::vector<std::pair<signed long, BaseSpinePtr>>> (GetNoSourceNeurons());
    DelayDDistribution = std::vector<std::vector<int>> (GetNoSourceNeurons());
    JCouplingDistribution = std::vector<std::vector<double>> (GetNoSourceNeurons());
    // hasPlasticityModel=targetPop->HasPlasticityModel();
    if (targetPopID == sourcePopID){
        isRecurrent=true;
    }
}

std::string Synapse::GetIdStr() const{
    NeuronInt sourcePopID  = this->sourcePop->GetId();
    NeuronInt targetPopID = this->targetPop->GetId();
    std::string idStr = std::to_string(targetPopID) + std::to_string(sourcePopID);
    return idStr;
}

std::string Synapse::GetIdStrWithULine() const {
    NeuronInt sourcePopID = this->sourcePop->GetId();
    NeuronInt targetPopID = this->targetPop->GetId();
    std::string idStr = std::to_string(targetPopID) + "_" + std::to_string(sourcePopID);
    return idStr;
}

void Synapse::LoadParameters(const std::vector<FileEntry>& synapseParameters){
    //std::cout << "loading parameters for synapse " << std::to_string(this->neuronsPre->GetId()) <<"->" <<std::to_string(this->neuronsPost->GetId()) << "\n";
    // These boolean variables keeps track of whether J_pot and P_pot are specified in the parameter file. If not, they are set to J and 0 respectively after the loop.
    bool foundJPot = false;

    //This should allow Synapse elements deletion
    this->isConnectedBool=true;

    for(auto& [parameterName, parameterValues] : synapseParameters) {
        if(parameterName.find("J_pot") != std::string::npos){
            this->Jpot = std::stod(parameterValues.at(0));
            foundJPot = true;
            //std::cout << "assigning J_pot = " << std::to_string(this->J_pot) << "\n";
        } else if(parameterName.find('J') != std::string::npos){
            this->J = std::stod(parameterValues.at(0));
            //std::cout << "assigning J = " << std::to_string(this->J) << "\n";
        } else if (parameterName.find("Sigma_j") != std::string::npos) {
            this->sigmaJ = std::stod(parameterValues.at(0));
            //std::cout << "assigning SigmaJ = " << std::to_string(this->SigmaJ) << "\n";
        } else if(parameterName.find("P_pot") != std::string::npos){
            this->Ppot = std::stod(parameterValues.at(0));
            //std::cout << "assigning P_pot = " << std::to_string(this->P_pot) << "\n";
        } else if(parameterName.find("D_min") != std::string::npos){
            this->Dmin = static_cast<int>(std::round(std::stod(parameterValues.at(0)) / infoGlobal->dtTimestep));
        } else if(parameterName.find("D_max") != std::string::npos){
            this->Dmax = static_cast<int>(std::round(std::stod(parameterValues.at(0)) / infoGlobal->dtTimestep));
        } else if(parameterName.find("connectivity_type") != std::string::npos){
            if(parameterValues.at(0) == IDstringRandomConnectivity || parameterValues.at(0) == IDstringHeteroRandomConnectivity){
                geometry = std::make_unique<RandomConnectivity>(this,infoGlobal);
            } else if(parameterValues.at(0) == IDstringPoissonConnectivity){
                geometry = std::make_unique<PoissonConnectivity>(this,infoGlobal);
            } else if (parameterValues.at(0) == IDstringDistanceConnectivity) {
				geometry = std::make_unique<DistanceConnectivity>(this, infoGlobal);
			// } else if(parameterValues.at(0) == IDstringHeteroRandomConnectivity) {
            //     this->geometry = std::make_unique<HeteroRandomConnectivity>(this, infoGlobal);
            } else if (parameterValues.at(0) == IDstringAdjacencyMatrixConnectivity) {
                this->geometry = std::make_unique<AdjacencyMatrixConnectivity>(this, infoGlobal);
            } else {
                throw "Connectivity type was not understood";
            }
        } else if(parameterName.find("connected") != std::string::npos){
            if (parameterValues.at(0).find("false")!=std::string::npos){
                this->isConnectedBool=false;
            } else if (parameterValues.at(0).find("true")!=std::string::npos){
                this->isConnectedBool=true; //just to make sure
            } else {
                throw "Unspecified bool for synapse_connected"; //if connected is there, but the input is not correct, exception
            }
        } else if(parameterName.find("targetBranch") != std::string::npos){
            if (parameterValues.at(0).find("random") != std::string::npos){
                this->branchTarget.randomTargetBranch=true;
            } else if (parameterValues.at(0).find("none") != std::string::npos){
                //Here nothing is done to handle the case where we do not used branched
                //Checking if the input is an integer would take more time than checking if None is in a single string with 5 characters max (I think).
            }else if (parameterValues.at(0).find("ordered") != std::string::npos){
                this->branchTarget.orderedTargetBranch=true;
            } else {
            this->branchTarget.targetBranch = std::stoi(parameterValues.at(0));
            if (this->branchTarget.targetBranch>=targetPop->GetNoBranches()){
                throw "Targeted branch does not exist";
            }
            this->branchTarget.setTargetBranch=true;
            //Missing exception management for when the input is not an integer.
            }
        } else if(parameterName.find("subregion") != std::string::npos){
            this->branchTarget.DendriticSubRegion = parameterValues.at(0).at(0); //Re-think this char array wacky stuff
        } else if(parameterName.find("relativeCoupling") != std::string::npos){
            this->relativeCouplingStrength = std::stod(parameterValues.at(0));
        } else if(parameterName.find("seed") != std::string::npos){
            SetSeed(std::stoi(parameterValues.at(0)));
            userSeed=true;
        }
    }

    // Set default parameterValues of the optional parameters J_pot and P_pot if necessary

    if (!foundJPot) {
		this->Jpot = this->J;
	}

	waitingMatrix.resize(GetNoTargetNeurons());
	for (NeuronInt targetNeuron : std::ranges::views::iota(0, GetNoTargetNeurons()))
		waitingMatrix.at(targetNeuron).resize(static_cast<size_t>(Dmax) + 1, 0.0);
    std::vector<FileEntry> connectivityParameters{FilterStringEntries(synapseParameters, "connectivity")};
    geometry->LoadParameters(connectivityParameters);
}

void Synapse::SaveParameters(std::ofstream& wParameterStream,std::string idString) const{
    wParameterStream << idString << "type\t\t\t\t\t\t" << GetTypeStr() << "\n";
    wParameterStream << idString << "connected\t\t\t\t\t\t" << std::boolalpha << this->isConnectedBool <<std::noboolalpha<< "\n";
    wParameterStream << idString << "D_min\t\t\t\t\t\t" << std::to_string(this->Dmin*infoGlobal->dtTimestep) << " #secs\n";
    wParameterStream << idString << "D_max\t\t\t\t\t\t" << std::to_string(this->Dmax*infoGlobal->dtTimestep) << " #secs\n";
    if (!this->ignoreJDistribution){
        wParameterStream << idString << "J\t\t\t\t\t\t\t" << std::to_string(this->J) << " #dmV/Spike\n";
        wParameterStream << idString << "Sigma_j\t\t\t\t\t\t" << std::to_string(this->sigmaJ) << " #dmV/Spike\n";
        wParameterStream << idString << "J_pot\t\t\t\t\t\t" << std::to_string(this->Jpot) << " #dmV/Spike\n";
        wParameterStream << idString << "P_pot\t\t\t\t\t\t" << std::to_string(this->Ppot) << "\n";
    }
    if (userSeed){
        wParameterStream << idString << "seed\t\t\t\t\t\t" << std::to_string(this->seed)  << "\n";
    }

    //Synaptic targeting
    if (this->targetPop->IsBranchedBool()){
        wParameterStream << idString << "targetBranch\t\t\t\t\t";
        if (this->branchTarget.randomTargetBranch){
            wParameterStream<<"random"; //Missing comments on what this is supposed to do
        } else if(this->branchTarget.setTargetBranch){
            wParameterStream<<std::to_string(this->branchTarget.targetBranch);//Missing comments on what this is supposed to do
        }else if(this->branchTarget.orderedTargetBranch){
        wParameterStream<<"ordered";
        } else {
            wParameterStream<<"none";//Missing comments on what this is supposed to do
        }
        wParameterStream << "\t\t\t#You can target branches in an 'ordered' manner (0,1,2...), 'random', or set (if you input a number). Put none if the HS does not used branched morphology\n";
    }
    if (this->ignoreJDistribution){
        wParameterStream << idString << "relativeCoupling\t\t\t\t" << std::to_string(relativeCouplingStrength)<< "\t#Relative coupling strength (only non-J plasticity models)";
    }
    if (geometry != nullptr){
        geometry->SaveParameters(wParameterStream,idString);
    }
}
void Synapse::AllocateSynapseStandard(NeuronInt targetNeuron, NeuronInt sourceNeuron) {
    targetSpineList.at(sourceNeuron).push_back(std::make_pair(targetNeuron, nullptr));
}

void Synapse::AllocateSynapseWithPlasticity(NeuronInt targetNeuron, NeuronInt sourceNeuron) {
    BaseSpinePtr synapseSpinePtr = targetPop->AllocateNewSynapse(targetNeuron, branchTarget);//everything except first var can be moved to syn ref

    if (synapseSpinePtr != nullptr) {
        if (ignoreJDistribution){
            synapseSpinePtr->SetRelativeCouplingStrength(relativeCouplingStrength);
        } else {
            //The SetRelativeCouplingStrength is already done in the PostConnectFunction (unscaled as of now), as J distribution is not written yet in flow
        }
        synapseSpinePtr->SetPreNeuronId(sourceNeuron);
        synapseSpinePtr->SetPostNeuronId(targetNeuron);

    } else {
        throw "AllocateNewSynapseWithPlasticity returned a nullptr";
    }
    targetSpineList.at(sourceNeuron).push_back(std::make_pair(targetNeuron, synapseSpinePtr)); //Const cast here too
    //synapticTargets.at(sourceNeuron).push_back(targetNeuron);//Are both truly necessary?
}

// double Synapse::GetRecurrentInput(NeuronInt targetNeuron) const{
// 	return waitingMatrix.at(targetNeuron).at((this->infoGlobal->timeStep) % (Dmax + 1));
// }

void Synapse::Advect(std::vector<double>& synaptic_dV, std::mutex& _syndVMutex) {
    //REWRITE this function. There is no need to resize and refill, simply input as rvalue in FillWaitingMatrix. This requires a re-write of all advectSpikers sadly
    ResetcumulatedDV();//Does this do something ?


    //Get list of spikers
    const std::vector<NeuronInt>& spikersUsed = sourcePop->GetSpikers();

    //Go through all the spikers and add current arising from spikers to waiting_matrix //PARALLELIZE EXECUTION HERE
    // std::for_each(spikersUsed.begin(), spikersUsed.end(), [this](NeuronInt spiker){
    //     FillWaitingMatrix(spiker, AdvectSpikers(spiker));//Here target list is used again OPTIMIZATION
    // });
    std::for_each(spikersUsed.begin(), spikersUsed.end(), [this](NeuronInt spiker){
        FillWaitingMatrix(spiker, AdvectSpikers(spiker));//Here target list is used again OPTIMIZATION
    });
    // for(NeuronInt spiker: spikersUsed){

    //     FillWaitingMatrix(spiker, AdvectSpikers(spiker));//Here target list is used again OPTIMIZATION
    // }

	ReadWaitingMatrixEntry(synaptic_dV,_syndVMutex);
    //std::cout << cumulatedDV << std::endl;
	//advect_finalize(&waiting_matrix);// Is it OK for Conductance, Exponential and Graupner Synapse ? Fill Waiting Matrix has already been called
}

void Synapse::FillWaitingMatrix(NeuronInt spiker, std::vector<double>&& currents){
    //std::vector<signed long> synapticTargets {geometry->GetTargetList(spiker)};
    //std::vector<int> *tD = geometry->GetDistributionD(spiker);
    std::vector<std::pair<NeuronInt, BaseSpinePtr>>& singleTargetList = targetSpineList.at(spiker);
    for(NeuronInt targetNeuronIndex : std::ranges::views::iota(0, static_cast<NeuronInt>(singleTargetList.size()))){
        // int delay {GetDistributionD(targetNeuronIndex,spiker)};                              // get individual synaptic delay between spiker and target
        size_t matrixIndex {(static_cast<size_t>(this->infoGlobal->timeStep) + static_cast<size_t>(GetDistributionD(targetNeuronIndex,spiker)))%static_cast<size_t>(Dmax+1)};
        //std::cout<<delay<<"-"<<matrixIndex<<"-"<<waitingMatrix.at(synapticTargetList.at(spiker).at(targetNeuronIndex).first).size()<<"-"<<targetNeuronIndex<<"\n";
        //std::lock_guard<std::mutex> _lockMutex(_waitingMatrixMutexLock);
        waitingMatrix.at(singleTargetList.at(targetNeuronIndex).first).at(matrixIndex) += currents.at(targetNeuronIndex); // add spiker to waiting matrix
    }
}

double Synapse::GetCouplingStrength(NeuronInt targetNeuronIndex, NeuronInt sourceNeuron) const {
    //Here goes everything inside GetCouplingStrength, but plasticity models will not work hand in hand with network scaling (local or global)
    if(infoGlobal->networkScaling_mode == 0){
        return GetDistributionJ(targetNeuronIndex, sourceNeuron)*localScalingFactor;
    } else if(infoGlobal->networkScaling_mode == 1){
        return GetDistributionJ(targetNeuronIndex, sourceNeuron)*globalScalingFactor;
    } else {
        throw "error in Synapse::PostCOnnectNeurons()";
    }
}

void Synapse::SetSeed(int inputSeed) {
    seed       = inputSeed;
    generator  = std::mt19937(inputSeed);
}

void Synapse::ResetWaitingMatrixEntry() {
    long index {(this->infoGlobal->timeStep) % (Dmax + 1)};
    for (std::vector<double>& targetNeuronWaitingVector : waitingMatrix){
        targetNeuronWaitingVector.at(index)=0;
    }
}

void Synapse::ReadWaitingMatrixEntry(std::vector<double>& synaptic_dV, std::mutex& _syndVMutex) const {
	long index {(this->infoGlobal->timeStep) % (Dmax + 1)};
    std::lock_guard<std::mutex> _lockMutex(_syndVMutex);
	for (NeuronInt targetNeuron : std::ranges::views::iota(0, GetNoTargetNeurons())){
		synaptic_dV.at(targetNeuron) += waitingMatrix.at(targetNeuron).at(index);//This is the only non-threadsafe line for par, while par_unseq I think it is best to leave alone
    }
}

void Synapse::ConnectNeurons(){
    if (IsConnected()){
        geometry->ConnectNeurons();
    } else {
        std::cout<< "Connection skipped." <<"\n";
    }
    localScalingFactor=pow(geometry->GetExpectedConnections(),infoGlobal->networkScaling_synStrength);
    globalScalingFactor=pow(infoGlobal->totalNeurons,infoGlobal->networkScaling_synStrength);
}

void Synapse::PostConnectNeurons() {
    for (std::vector<std::pair<NeuronInt, BaseSpinePtr>>& singleSourceList : targetSpineList){
        //Here we sort the subentries of targetlist so that multiple synapses do not generate conflicts
        std::sort(singleSourceList.begin(), singleSourceList.end(), [](std::pair<NeuronInt, BaseSpinePtr> pair1,std::pair<NeuronInt, BaseSpinePtr> pair2){return pair1.first<pair2.first;});
    }
    // if (HasJDistribution){
    //     if(infoGlobal->networkScaling_mode == 0){
    //         for (NeuronInt sourceNeuron : std::ranges::views::iota(0,sourcePop->GetNoNeurons())){
    //             for (NeuronInt targetNeuron : std::ranges::views::iota(0,targetPop->GetNoNeurons())){
    //                 JCouplingDistribution.at(sourceNeuron).at(targetNeuron)*=localScalingFactor;
    //             }
    //         }
    //     }
    //     else if(infoGlobal->networkScaling_mode == 1){
    //         for (NeuronInt sourceNeuron : std::ranges::views::iota(0,sourcePop->GetNoNeurons())){
    //             for (NeuronInt targetNeuron : std::ranges::views::iota(0,targetPop->GetNoNeurons())){
    //                 JCouplingDistribution.at(sourceNeuron).at(targetNeuron)*=globalScalingFactor;
    //             }
    //         }
    //      } else {
    //         throw "error in Synapse::PostCOnnectNeurons()";
    //         }
    // } else {
    //     if(infoGlobal->networkScaling_mode == 0){
    //         scaledJ=J*localScalingFactor;
    //     } else if(infoGlobal->networkScaling_mode == 1){
    //         scaledJ=J*globalScalingFactor;
    //     } else {
    //         throw "error in Synapse::PostCOnnectNeurons()";
    //     }
    // }
    if ((!hasPlasticityModel && targetPop->HasPlasticityModel()) || (hasPlasticityModel && !targetPop->HasPlasticityModel())){
        throw "Logical error in the usage of plasticity models.";
    }
    if (!ignoreJDistribution && hasPlasticityModel){
        for(NeuronInt sourceNeuron : std::ranges::views::iota(0, static_cast<NeuronInt>(targetSpineList.size()))){
            for(auto& [targetNeuron, spinePtr] : targetSpineList.at(sourceNeuron)){
                spinePtr->SetRelativeCouplingStrength(GetCouplingStrength(targetNeuron,sourceNeuron)); //const cast here
            }
        }
    }
    //Here put all scaling things(or before this line) to avoid calculating them in GetCouplingStrength
    this->targetPop->PostConnectSetUp();//For copying of pointers to synaptic branches
}

//Connectivity functions
void Synapse::WriteConnectivity(const std::string &filename, NeuronInt noRecordedNeuronsConnectivity) const {

    if(noRecordedNeuronsConnectivity <= 0){
        return;
    }

    std::ofstream connectivityFileStream(filename + ".txt");
    
    int count;
    size_t writtenIndex;
    NeuronInt noTargetNeurons = std::min(noRecordedNeuronsConnectivity, GetNoTargetNeurons());
    // NeuronInt noSourceNeurons = std::min(noRecordedNeuronsConnectivity, GetNoSourceNeurons());

    //iterate through sourceNeuron & target neurons & write connectivity
    for (const std::vector<std::pair<NeuronInt, BaseSpinePtr>>& singleSourceList : targetSpineList){
        writtenIndex = 0;
        for(NeuronInt targetNeuron : std::ranges::views::iota(0, noTargetNeurons)){
			if (singleSourceList.size() <= writtenIndex || singleSourceList.at(writtenIndex).first != targetNeuron){
				connectivityFileStream << "0 ";
			} else {
				count = 1;
				writtenIndex++;
				while (singleSourceList.size() > writtenIndex && singleSourceList.at(writtenIndex).first == targetNeuron) { //in case the same pair can be connected multiple times
					writtenIndex++;
					count++;
				}
				connectivityFileStream << std::to_string(count)<<" ";
			}
        }
        connectivityFileStream << "\n";

        /*stream << "Targets for neuron " << sourceNeuron << ":" << synapticTargets.at(sourceNeuron).size() << " -- ";
        for(auto const& target: (synapticTargets.at(sourceNeuron))){
            stream << std::to_string(target) << " ";
        }
        stream << "\n";*/
    }
    connectivityFileStream.close();

}
void Synapse::SetDistributionD() {
    std::uniform_int_distribution<int> uniformDistribution(this->GetMinD(),this->GetMaxD());
	NeuronInt noTargets;
	if (GetMaxD()== GetMinD()) {
		HasDdistribution = false;
		return;
	}
    //std::cout << "printing D_distribution: \n";
    for(NeuronInt sourceNeuron : std::ranges::views::iota(0,GetNoSourceNeurons())){
        if(!targetSpineList.at(sourceNeuron).empty()){
            noTargets = static_cast<NeuronInt>(GetNoTargetedNeurons(sourceNeuron));
			DelayDDistribution.at(sourceNeuron).resize(noTargets, 0);
			for (NeuronInt targetNeuron : std::ranges::views::iota(0,noTargets)) {
                DelayDDistribution.at(sourceNeuron).at(targetNeuron)=uniformDistribution(generator);
                //std::cout << std::to_string(d) << " ";
            }
            //std::cout << "\n";
        }
    }
}
void Synapse::SetDistributionJ() {
    // double rdj;//random number to determine if the synapse is potentiated
    // double rds;//random number to determine the deviation to the mean j

	if (this->GetSigmaJ() <= std::numeric_limits<double>::epsilon() && (this->GetPPot() <= std::numeric_limits<double>::epsilon() || this->GetJPot() == this->GetJBase())) {
        HasJDistribution = false;
		return;
	}
    std::uniform_real_distribution<double> realDistribution(0,1);
    std::normal_distribution<double> gaussianDistribution(0, 1);
	NeuronInt noTargets;
    double JCouplingStrength;
    //std::cout << "printing J_distribution: \n";
    for(NeuronInt sourceNeuron : std::ranges::views::iota(0, this->GetNoSourceNeurons())){
        noTargets = static_cast<NeuronInt>((targetSpineList.at(sourceNeuron)).size());
		JCouplingDistribution.at(sourceNeuron).resize(noTargets, 0.0);
		//std::cout <<  "\n";
        for(NeuronInt targetNeuron : std::ranges::views::iota(0, noTargets)){
			// rdj = real_distribution(generator);
            // rds = Gaussian(generator);
            JCouplingStrength = (realDistribution(generator) < this->GetPPot())?this->GetJPot():this->GetJBase();
            JCouplingStrength += this->GetSigmaJ() * gaussianDistribution(generator);
            //if ((this->GetJPot() * this->GetJBase() >= 0) && (j * this->GetJBase() < 0)) {
            //    j = 0;//Dale's law will be strictly enforced unless the provided J and Jpot already blatantly violate iterator
            //}
            JCouplingDistribution.at(sourceNeuron).at(targetNeuron) = JCouplingStrength;
		}
		//std::cout << "\n";
    }
}
void Synapse::WriteDistributionD(const std::string& filename, NeuronInt noRecordedNeuronsDelay) const {
    if(noRecordedNeuronsDelay <= 0){
        return;
    }
    int count;
    size_t writtenConnection;
    double delay;
    std::ofstream DDistributionFileStream(filename + ".txt"); //ios::binary

    //iterate through sourceNeuron & target neurons & write connectivity
    for(NeuronInt sourceNeuron : std::ranges::views::iota(0, std::min(GetNoSourceNeurons(), noRecordedNeuronsDelay))){
        writtenConnection = 0; // counts the number of connections that have been written for this sourceNeuron neuron
        for(NeuronInt targetNeuron : std::ranges::views::iota(0, std::min(GetNoTargetNeurons(), noRecordedNeuronsDelay))){
			if (targetSpineList.at(sourceNeuron).size() <= writtenConnection || targetSpineList.at(sourceNeuron).at(writtenConnection).first != targetNeuron){ // if all connections have been written OR target id does not match
				DDistributionFileStream << "nan ";
			} else {
				count = 1;
                delay = (HasDdistribution)? DelayDDistribution.at(sourceNeuron).at(writtenConnection) : GetMaxD();
				writtenConnection++;
				while (targetSpineList.at(sourceNeuron).size() > writtenConnection && targetSpineList.at(sourceNeuron).at(writtenConnection).first == targetNeuron) { // while connections are still to be written and target id matches
					writtenConnection++;
					count++;
                    delay += (HasDdistribution)? DelayDDistribution.at(sourceNeuron).at(writtenConnection) : GetMaxD();
				}
				DDistributionFileStream << std::to_string(delay/count)<<" ";
			}
        }
        DDistributionFileStream << "\n";

        /*stream << "Targets for neuron " << sourceNeuron << ":" << synapticTargets.at(sourceNeuron).size() << " -- ";
        for(auto const& target: (synapticTargets.at(sourceNeuron))){
            stream << std::to_string(target) << " ";
        }
        stream << "\n";*/
    }
    DDistributionFileStream.close();
}
void Synapse::WriteDistributionJ(const std::string& filename,NeuronInt noRecordedNeuronsJPot) const {
    if(noRecordedNeuronsJPot <= 0){
        return;
    }

    //This is a recorder-called function (after jumping through interfaces). The number of recorded neurons dictates what is written.
    //int count;
    size_t writtenConnection;
    double JCouplingStrength;
    std::ofstream JDistributionFileStream(filename + ".txt"); //ios::binary
    //iterate through source & target neurons & write connectivity
    for(NeuronInt sourceNeuron : std::ranges::views::iota(0, std::min(GetNoSourceNeurons(),noRecordedNeuronsJPot))){
        writtenConnection = 0; // counts the number of connections that have been written for this source neuron
        for(NeuronInt targetNeuron : std::ranges::views::iota(0,std::min(GetNoTargetNeurons(),noRecordedNeuronsJPot))){
			if (targetSpineList.at(sourceNeuron).size() <= writtenConnection || targetSpineList.at(sourceNeuron).at(writtenConnection).first != targetNeuron){ // if all connections have been written OR target id does not match
				//Here unless we finish (in which case we should not write anything) or there is a logical error (which seems likely from implementation) this will not run
                JDistributionFileStream << "nan ";
			} else {
                //If it all makes sense, we get J and set count to 1
				//count = 1;
                JCouplingStrength = (HasJDistribution)? JCouplingDistribution.at(sourceNeuron).at(writtenConnection) : GetJBase();
                //Now we increase wC count bc we are going to at least write J/1 in the file.
				writtenConnection++;
				while (targetSpineList.at(sourceNeuron).size() > writtenConnection && targetSpineList.at(sourceNeuron).at(writtenConnection).first == targetNeuron) { // while connections are still to be written and target id matches
					//And write every successive synapse between source and target (which will exclude non-successive, but it is fixable by sorting target list by first element.)
                    writtenConnection++;
					//count++;
                    JCouplingStrength = (HasJDistribution)? JCouplingDistribution.at(sourceNeuron).at(writtenConnection) : GetJBase();
				}
				JDistributionFileStream << std::to_string(JCouplingStrength)<<" ";
                //Why is this averaged?? Average coupling strength?
			}
        }
        JDistributionFileStream << "\n";

        /*stream << "Targets for neuron " << sourceNeuron << ":" << synapticTargets.at(sourceNeuron).size() << " -- ";
        for(auto const& target: (synapticTargets.at(sourceNeuron))){
            stream << std::to_string(target) << " ";
        }
        stream << "\n";*/
    }
    JDistributionFileStream.close();

}

void Synapse::AllocateSynapse(NeuronInt targetNeuron, NeuronInt sourceNeuron) {
    AllocateSynapseStandard(targetNeuron,sourceNeuron);
}
