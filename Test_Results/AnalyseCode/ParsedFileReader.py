import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

class NeuronPopRecording: #Suppose to be one file
    metadata_datanames=['noNeurons', 'dt', 'totalTimesteps', 'neuronPopId', 'simTime']
    def __init__(self, neuronPop_id : int):
        self.neuronPop_id=neuronPop_id
    def append_metadata(self, metadata_list : list):
        self.metadata=metadata_list
    def append_spiketimes(self, N_index : int, spiketimes_list : list):
        setattr(self, "Neuron_%s" % (N_index), spiketimes_list)
    def create_firing_matrix(self):
        self.F_matrix=np.zeros(shape=(self.metadata[2], self.metadata[0]))
        for neuron in dir(self):
            if neuron.startswith("Neuron"):
                neuron_index=int(neuron.split('_')[1])
                spiketimes=self.__dict__[neuron]
                for st in spiketimes:
                    self.F_matrix[int(st/self.metadata[1])-1, neuron_index]=1

class OutputRecordingList: #Compendium of all recordings
#    def __iter__(self):
#        for attr, value in self.__dict__.iteritems():
#            yield attr, value
    def get(self):
        return self.__dict__
    def append(self, var : NeuronPopRecording):
        setattr(self, "NeuronPop_%s" % (var.neuronPop_id), var)

#Here we use split = to get the name and the values, then split the second with commas to get the spike-time list
def read_reasterplot_file_into_object(neuronPop_files):
    all_recordingsST=OutputRecordingList()
    for i, file in enumerate(neuronPop_files):
        with open(file, "r") as file:
            lines=file.readlines()
            for line in lines:
                if line.startswith('M'):
                    subline=line.split('=')
                    neuronpop_id=int(subline[0].split('_')[1])
                    metadata_list=subline[1].split(',')

                    metadata_list[0]=int(metadata_list[0])
                    metadata_list[1]=float(metadata_list[1])
                    metadata_list[2]=int(metadata_list[2])
                    metadata_list[3]=int(metadata_list[3])
                    metadata_list[4]=float(metadata_list[4])

                    recording=NeuronPopRecording(neuronpop_id)
                    recording.append_metadata(metadata_list)
                    all_recordingsST.append(recording)

                elif line.startswith('N'):
                    subline=line.split('=')
                    idline=subline[0].split('_')
                    neuronpop_id=int(idline[1])
                    neuron_id=int(idline[2])

                    spiketimes=list()
                    for i in subline[1].split(','):
                        if i=='':
                            continue
                        elif i=='\n':
                            continue
                        else:
                            spiketimes.append(float(i))
                    all_recordingsST.__getattribute__(list(all_recordingsST.__dict__.keys())[neuronpop_id]).__getattribute__("append_spiketimes").__call__(neuron_id, spiketimes)
                else:
                    print('Oversight error')
                    raise Exception
        all_recordingsST.append(recording)
    return all_recordingsST

### Parsed spike times file reading

#File reading here (files are always iterated on alphabetical order)
parsedOutput_string='_ParsedOutput.dat'

base_dir = os.getcwd()

neuronPop_files = [f for f in os.listdir(base_dir) if f.endswith(parsedOutput_string)]

### Reading parsed spike times file (Neuron output) and creating firing matrix
all_recordings=read_reasterplot_file_into_object(neuronPop_files)

for attr in all_recordings.__dict__:
    print(attr)
    all_recordings.__getattribute__(attr).__getattribute__("create_firing_matrix").__call__()
    ##Up until here everything woks as expected

### Rasterplot parsing code to spike times
# (They will be lacking non-firing neurons, and there is no metadata)
#File reading here (files are always iterated on alphabetical order)
rasterplot_string='Rasterplot.dat'

base_dir = os.getcwd()

rasterplot_files = [f for f in os.listdir(base_dir) if f.endswith(rasterplot_string)]
rasterdata=pd.read_table(rasterplot_files[0], comment='#')
neuronPopIds=set(rasterdata['Pop_id'])
neuronIds=set(rasterdata['Neuron_id'])
all_recordings2=OutputRecordingList()
for i in neuronPopIds:
    all_recordings2.append(NeuronPopRecording(i))
    for j in neuronIds:
        spiketimes = rasterdata[rasterdata['Neuron_id']==j]
        spiketimes= spiketimes[spiketimes['Pop_id']==i]
        spiketimes = list(spiketimes["Spike_t"])
        all_recordings2.__getattribute__(list(all_recordings2.__dict__.keys())[i]).__getattribute__("append_spiketimes").__call__(j, spiketimes)
