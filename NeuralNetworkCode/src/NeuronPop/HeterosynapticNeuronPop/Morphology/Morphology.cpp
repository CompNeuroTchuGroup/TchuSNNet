#include "Morphology.hpp"

Morphology::Morphology(GlobalSimInfo *info): info(info), lastPostSpikeTime(-200), weightsSum(0), totalPostSpikes(0), totalPreSpikes(0), weightNormalization(NOPNormalization) {
}

void Morphology::LoadParameters(std::vector<std::string>* input) {

    std::string name;
    std::vector<std::string> values;

    // checks for correct initialization

    for (auto& it : *input) {
        SplitString(&it, &name, &values);

        if (name.find("seed") != std::string::npos) {
            this->seed = std::stoi(values.at(0));
            this->generator=std::mt19937(this->seed);
        }
        //include here max and min weights
        if(info->globalSeed != -1){
        std::uniform_int_distribution<int> distribution(0,INT32_MAX);
        this->seed = distribution(info->globalGenerator);
        this->generator = std::mt19937(seed);
    }
    }
}
void Morphology::SaveParameters(std::ofstream *stream, std::string neuronPreId) {
    *stream<< "#########From here on is all Heterosynaptic plasticity###########\n";

    *stream << neuronPreId<<"_morphology_type\t\t\t"<<this->GetType()<<"\n";

    *stream << neuronPreId<<"_morphology_seed\t\t\t"<<this->seed<<"\n";
}

void Morphology::RecordPostSpike() {
    this->lastPostSpikeTime = this->info->dt * static_cast<double> (this->info->time_step);
    this->totalPostSpikes++;
    this->postSpiked = true;
}


void Morphology::RecordExcitatoryPreSpike(int spikedSynapseId) {
        //This function is NOT DELAY COMPATIBLE (careful with the delays in synapse objects)
    //Is there supposed to be a different Inhibitory function?
    this->totalPreSpikes++;
    // STDP Analysis
    //this->preSpikes.emplace_back(spikedSynapseId, this->baseSynapseData.at(spikedSynapseId)->lastSpike);
}

std::valarray<double> Morphology::GetIndividualSynapticProfile(unsigned long synapseId) const
{
    return baseSynapseData.at(synapseId)->GetIndividualSynapticProfile();
}

std::valarray<double> Morphology::GetOverallSynapticProfile()
{
    /*
     * returned array organised as follows:
     * item 1: average synaptic weight
     * item 2: total post spikes
     * item 3: total pre spikes
     * item 4: average plasticity events
     * */
    std::valarray<double> dataArray(4);
    size_t sizeOfSynapseData {this->baseSynapseData.size()};
    double weightSum = std::accumulate(this->baseSynapseData.begin(), this->baseSynapseData.end(), 0.0,
                                       [] (double acc, const BaseSpinePtr& syn) { return acc + syn->GetWeight(); });
    //CalcMorphoPlasticityEvents();

   dataArray[0] = weightSum / sizeOfSynapseData;
   dataArray[1] = this->totalPostSpikes;
   dataArray[2] = this->totalPreSpikes;
   dataArray[3] = static_cast<double>(totalPlasticityEvents) / sizeOfSynapseData;
    return dataArray;
}

std::string Morphology::GetOverallSynapticProfileHeaderInfo() const
{
    return std::string("{<average weight>, <total post spikes>, <total pre spikes>, <average number of plasticity events>}");
}

std::string Morphology::GetIndividualSynapticProfileHeaderInfo() const
{
    return baseSynapseData.at(0)->GetIndividualSynapticProfileHeaderInfo();
}

unsigned long Morphology::GetSynapseCount() const {
    return static_cast<unsigned long>(this->baseSynapseData.size());
}

double Morphology::GetWeight(unsigned long synapseId) const {
    return this->baseSynapseData.at(synapseId)->GetWeight();
}

void Morphology::NormalizeWeights() {
    if (this->weightNormalization == HardNormalization) {
        this->hardNormalize();
    } else if (this->weightNormalization == SoftMaxNormalization) {
        this->softMaxNormalize();
    }
}

void Morphology::hardNormalize() {
    for (auto& syn: this->baseSynapseData) {
        syn->SetWeight(std::max(minWeight, std::min(maxWeight, syn->GetWeight())));
    }
}

void Morphology::softMaxNormalize() {
    double maxWeights, sumWeights;

    maxWeights = std::numeric_limits<double>::min();
    for (auto& syn : this->baseSynapseData) {
        maxWeights = std::max(maxWeights, syn->GetWeight());
    }

    sumWeights = 0;
    for (auto& syn : this->baseSynapseData) {
        sumWeights += std::exp(syn->GetWeight() - maxWeights);
    }
}

void Morphology::WeightDecay() {
    if (this->decayWeights) {
        for (const BaseSpinePtr& syn: this->baseSynapseData) {
            syn->SetWeight(syn->GetWeight() * weightExpDecay);
        }
    }
}

double Morphology::GenerateSynapticWeight(){
    double weight{};
    std::uniform_real_distribution<double> distribution(this->minWeight,this->maxWeight);
            if (this->distributeWeights) {
                std::mt19937& generatorRef = this->generator;
                weight = distribution(generatorRef);
            } else {
                weight = this->initialWeights; // assuming a range of weight between 0 and 2, weight is initialized to midpoint: 1
            }
        //this->weightsSum += weight;
        return weight;
}
/*
void Morphology::triggerStatOut(std::string dirPath) {
    std::ofstream preFile;
    std::ofstream postFile;

    std::ofstream thetasFile;
    std::ofstream weightsFile;

    preFile.open (dirPath + "_pres.dat", std::ofstream::out | std::ofstream::trunc);
    postFile.open (dirPath + "_posts.dat", std::ofstream::out | std::ofstream::trunc);
    thetasFile.open (dirPath + "_thetas.dat", std::ofstream::out | std::ofstream::trunc);
    weightsFile.open (dirPath + "_weights.dat", std::ofstream::out | std::ofstream::trunc);

    for (const auto& line: preSpikes) {
        preFile << line.first << ":" << line.second << std::endl;
    }

    for (const auto& line: postSpikes) {
        postFile << line << std::endl;
    }

    unsigned long count = 0;
    for (const auto& line: theta_changes) {
        if (count % 2 == 0) {
            thetasFile << "s -- ";
        }
        thetasFile << line.first << ":" << line.second << std::endl;
        count++;
    }

    for (const auto& line: weight_changes) {
        weightsFile << line.first << ":" << line.second << std::endl;
    }

    preFile.close();
    postFile.close();
    thetasFile.close();
    weightsFile.close();

}*/

// void Morphology::printThetasAndWeights() {
//     std::cout << "weights: " << std::endl;
//     for (auto& syn: this->baseSynapseData) {
//         std::cout << syn->GetWeight() << ", ";
//     }
//     std::cout << std::endl;

//     std::cout << "thetas: " << std::endl;
//     for (auto& syn: this->baseSynapseData) {
//         std::cout << syn->GetTheta() << ", ";
//     }
//     std::cout << std::endl;
// }


