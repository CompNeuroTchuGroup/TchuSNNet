import numpy as np
import copy
import random as rd
import os


class InputGenerator:
    def __init__(
        self, N_neurons, total_time_intervals,total_time, testname="sim", pop_id1="0",first_interval_true=True,max_offset_inpop=0
    ) -> None:
        self.total_time_intervals=total_time_intervals
        self.interval_duration=total_time/total_time_intervals
        self.max_offset_inpop=max_offset_inpop
        self.filename = testname + "_DictatNeuronPop_" + pop_id1 + "_spikers.txt"
        self.N_neurons = N_neurons
        self.max_frequencies = np.zeros(N_neurons)
        self.neuron_offsets = np.zeros(N_neurons)
        # self.neuron_ids = np.array(range(0, N_neurons))
        self.instruction_start_times = np.zeros((N_neurons, total_time_intervals))
        self.instruction_end_times = np.zeros((N_neurons, total_time_intervals))
        self.frequencies = np.zeros((N_neurons, total_time_intervals))

    def _generate_random_arrays(self, neuron_id: int, max_frequency):
        ###This actually generates a time interval that determines the phase
        rd.seed(neuron_id)
        self.neuron_offsets[neuron_id] = rd.uniform(-self.max_offset_inpop/2, +self.max_offset_inpop/2)
        self.max_frequencies[neuron_id] = rd.uniform(max_frequency / 2, max_frequency)

    def _generate_nonrandom_arrays(self, neuron_id: int, max_frequency):
        ###This actually generates a time interval that determines the phase
        self.neuron_offsets[neuron_id] = 0
        self.max_frequencies[neuron_id] = max_frequency

    def set_up_arrays_random_max_frequency(
        self, max_frequency: float, randomized: bool
    ):
        if randomized:
            for i in range(0, self.N_neurons):
                self._generate_random_arrays(i, max_frequency)
        else:
            for i in range(0, self.N_neurons):
                self._generate_nonrandom_arrays(i, max_frequency)
    def _max_step(self,neuron,iteration):
        self.instruction_start_times[neuron,iteration]=iteration*self.interval_duration+self.neuron_offsets[neuron]
        self.instruction_end_times[neuron,iteration]=(iteration+1)*self.interval_duration+self.neuron_offsets[neuron]
        self.frequencies[neuron,iteration]=self.max_frequencies[neuron]
    def _zero_step(self,neuron,iteration):
        self.instruction_start_times[neuron,iteration]=iteration*self.interval_duration+self.neuron_offsets[neuron]
        self.instruction_end_times[neuron,iteration]=(iteration+1)*self.interval_duration+self.neuron_offsets[neuron]
        self.frequencies[neuron,iteration]=0
    def _alternate_step_zero_to_max(self,start:int,neuron):
        for j in range (0,self.total_time_intervals):
            if ((j+start)%2==0):
                self._max_step(neuron,j)
            else:
                self._zero_step(neuron,j)
    def zero_to_max_alterate_step_instructions_shuffled(self):
        for i in range (0,self.N_neurons):
            self._alternate_step_zero_to_max(i%2,i)
    def zero_to_max_alterate_step_instructions_non_shuffled(self):
        for i in range (0,self.N_neurons):
            self._alternate_step_zero_to_max(0,i)
    def _one_step(self,neuron:int,iteration:int, frequency:float):
        self.instruction_start_times[neuron,iteration]=iteration*self.interval_duration+self.neuron_offsets[neuron]
        self.instruction_end_times[neuron,iteration]=(iteration+1)*self.interval_duration+self.neuron_offsets[neuron]
        self.frequencies[neuron,iteration]=frequency
    def _two_step(self,neuron:int,iteration:int, frequency:float):
        self.instruction_start_times[neuron,iteration]=iteration*self.interval_duration+self.neuron_offsets[neuron]
        self.instruction_end_times[neuron,iteration]=(iteration+1)*self.interval_duration+self.neuron_offsets[neuron]
        self.frequencies[neuron,iteration]=frequency
    def _alternate_step_one_to_two(self,start:int,neuron, frequency1, frequency2):
        for j in range (0,self.total_time_intervals):
            if ((j+start)%2==0):
                self._one_step(neuron,j, frequency1)
            else:
                self._two_step(neuron,j, frequency2)
    def one_to_two_alterate_step_instructions_shuffled(self,frequency1:float, frequency2:float):
        for i in range (0,self.N_neurons):
            self._alternate_step_one_to_two(i%2,i, frequency1, frequency2)
    def one_to_two_alterate_step_instructions_non_shuffled(self,frequency1:float, frequency2:float):
        for i in range (0,self.N_neurons):
            self._alternate_step_one_to_two(0,i, frequency1, frequency2)
    def generate_instructions_from_arrays(self):
        with open (self.filename, 'w') as text_file:
            for i in range (0,self.N_neurons):
                for j in range (0,self.total_time_intervals):
                    text_file.write('> '+ str(i)+' '+str(self.instruction_start_times[i,j])+' '+str(self.instruction_end_times[i,j])+ ' ' + str(self.frequencies[i,j]))
                    if i+1 == self.N_neurons and j+1 == self.total_time_intervals:
                        pass
                    else:
                        text_file.write('\n')
        text_file.close()

class SimpleGenerator:
    def __init__(self, N_neurons,total_time, total_intervals, first_interval_true=False,testname="sim", pop_id1="0",max_offset_inpop=0) -> None:
        # self.interval_duration=total_time/total_intervals
        self.end_time=total_time
        self.intervals=total_intervals
        self.first=first_interval_true
        self.max_offset_inpop=max_offset_inpop
        self.filename  = f"{testname}/{testname}_DictatNeuronPop_{pop_id1}_spikers.txt"
        self.N_neurons = N_neurons
        self.max_frequencies = np.zeros(N_neurons)
        self.neuron_offsets = np.zeros(N_neurons)
        self.neuron_ids = np.array(range(0, N_neurons))
        self.instruction_start_times = np.zeros((N_neurons, total_intervals))
        self.instruction_end_times = np.zeros((N_neurons, total_intervals))
        self.frequencies = np.zeros((N_neurons, total_intervals))
    def _set_times(self,startTimes, endTimes):
        if self.intervals==len(startTimes) and len(endTimes)==len(startTimes):
            for i in range (0,self.N_neurons):
                self.instruction_start_times[i]=copy.copy(startTimes)
                self.instruction_end_times[i]=copy.copy(endTimes)
        else:
            raise Exception
    def _set_freqs_on_off(self,frequency):
        mod=0
        if not self.first:
            mod=1
        for i in range (0,self.N_neurons):
            for j in range (0,self.intervals):
              if (mod+j)%2==0:
                self.frequencies[i][j]=frequency
    def generate_instructions_from_arrays(self):
        if not os.path.exists(os.path.dirname(self.filename)):
            os.makedirs(os.path.dirname(self.filename))
        with open (self.filename, 'w') as text_file:
            for i in range (0,self.N_neurons):
                for j in range (0,self.intervals):
                    text_file.write('> '+ str(i)+' '+str(self.instruction_start_times[i,j])+' '+str(self.instruction_end_times[i,j])+ ' ' + str(self.frequencies[i,j]))
                    if i+1 == self.N_neurons and j+1 == self.total_time_intervals:
                        pass
                    else:
                        text_file.write('\n')
        text_file.close()
    def SW_stim_protocol_one_burst(self,stimStart, noSpikes, frequency,offset=0):
        if self.intervals!=3:
            raise Exception
        stimEnd=noSpikes*(1/frequency)
        startTimes=[0,stimStart+offset,stimEnd+offset]
        endTimes=[stimStart+offset,stimEnd+offset,self.end_time]
        self._set_times(startTimes,endTimes)
        self._set_freqs_on_off(frequency)
        self.generate_instructions_from_arrays()

    
        
class PairedPopsGenerator:
    def __init__(self, testname="testCorr2", pop_id1="0", pop_id2="1",max_offset_inpop=0,first_interval_true=False, N_Neurons_1=0,N_Neurons_2=0,total_time=0,time_intervals=0) -> None:
        self.gen1=SimpleGenerator(N_neurons=N_Neurons_1,testname=testname,pop_id1=pop_id1,first_interval_true=first_interval_true,max_offset_inpop=max_offset_inpop,total_time=total_time,total_intervals=time_intervals)
        self.gen2=SimpleGenerator(testname=testname,pop_id1=pop_id2,first_interval_true=first_interval_true,max_offset_inpop=max_offset_inpop,total_time=total_time,total_intervals=time_intervals,N_neurons=N_Neurons_2)
    def paired_stimm(self,stimStart,noSpikes, deltaT, frequency):
        if deltaT>=0:
          self.gen1.SW_stim_protocol_one_burst(stimStart,noSpikes,frequency,offset=deltaT)
          self.gen2.SW_stim_protocol_one_burst(stimStart,noSpikes,frequency)
        else:
          deltaT=-deltaT
          self.gen1.SW_stim_protocol_one_burst(stimStart,noSpikes,frequency)
          self.gen2.SW_stim_protocol_one_burst(stimStart,noSpikes,frequency,offset=deltaT)
            

#This functions are useless, they give spiketimes
# def allocate_pos_deltat(
#     timesteps, offset, steps, presynaptic, postsynaptic, deltaT_steps, noSpikes
# ):
#     counter = 0
#     counter2 = 0
#     flag = False
#     for i in timesteps:
#         if i < offset or noSpikes < 1:
#             continue
#         elif counter % steps == 0 and not flag:
#             presynaptic.append(i)
#             counter2 = 0
#             flag = True

#         elif counter2 % deltaT_steps == 0 and flag:
#             postsynaptic.append(i)
#             flag = False
#             noSpikes -= 1
#         counter2 += 1
#         counter += 1


# def allocate_neg_deltat(
#     timesteps, offset, steps, presynaptic, postsynaptic, deltaT_steps, noSpikes
# ):
#     counter = 0
#     counter2 = 0
#     flag = False
#     for i in timesteps:
#         if i < offset or noSpikes < 1:
#             continue
#         elif counter % steps == 0 and not flag:
#             postsynaptic.append(i)
#             counter2 = 0
#             flag = True

#         elif counter2 % deltaT_steps == 0 and flag:
#             presynaptic.append(i)
#             flag = False
#             noSpikes -= 1
#         counter2 += 1
#         counter += 1


# def allocate_zero_dt(
#     timesteps, offset, steps, presynaptic, postsynaptic, deltaT_steps, noSpikes
# ):
#     counter = 0
#     for i in timesteps:
#         if i < offset or noSpikes < 1:
#             continue
#         elif counter % steps == 0:
#             postsynaptic.append(i)
#             presynaptic.append(i)
#             noSpikes -= 1

#         counter += 1
