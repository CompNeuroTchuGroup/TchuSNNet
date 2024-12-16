import numpy as np
import random as rd

### Base parameters and arrays

testname='testCorr2'
pop_id='0'
filename=testname+'_DictatNeuronPop_'+pop_id+'_spikers.txt'
N_neurons=200
total_time=600.0 #In seconds
total_time_intervals=2
interval_duration=total_time/total_time_intervals

max_frequencies=np.zeros(N_neurons)
neuron_offsets=np.zeros(N_neurons)

neuron_ids=np.array(range(0,N_neurons))
instruction_start_times=np.zeros((N_neurons,total_time_intervals))
instruction_end_times=np.zeros((N_neurons,total_time_intervals))
frequencies=np.zeros((N_neurons,total_time_intervals))

### Array functions

def _generate_random_arrays(neuron_id:int, max_frequency):
    ###This actually generates a time interval that determines the phase
    rd.seed(neuron_id)
    neuron_offsets[neuron_id]=rd.uniform(0.0, 0.05)
    max_frequencies[neuron_id]=rd.uniform(max_frequency/2, max_frequency)
def _generate_nonrandom_arrays(neuron_id:int, max_frequency):
    ###This actually generates a time interval that determines the phase
    neuron_offsets[neuron_id]=0
    max_frequencies[neuron_id]=max_frequency

def set_up_arrays_random_max_frequency(max_frequency:float, randomized:bool):
    if (randomized):
        for i in range (0,N_neurons):
            _generate_random_arrays(i, max_frequency)
    else:
        for i in range (0,N_neurons):
            _generate_nonrandom_arrays(i, max_frequency)

### Instruction generation

def _max_step(neuron,iteration):
    instruction_start_times[neuron,iteration]=iteration*interval_duration+neuron_offsets[neuron]
    instruction_end_times[neuron,iteration]=(iteration+1)*interval_duration+neuron_offsets[neuron]
    frequencies[neuron,iteration]=max_frequencies[neuron]
def _zero_step(neuron,iteration):
    instruction_start_times[neuron,iteration]=iteration*interval_duration+neuron_offsets[neuron]
    instruction_end_times[neuron,iteration]=(iteration+1)*interval_duration+neuron_offsets[neuron]
    frequencies[neuron,iteration]=0
def _alternate_step_zero_to_max(start:int,neuron):
    for j in range (0,total_time_intervals):
        if ((j+start)%2==0):
            _max_step(neuron,j)
        else:
            _zero_step(neuron,j)
def zero_to_max_alterate_step_instructions_shuffled():
    for i in range (0,N_neurons):
        _alternate_step_zero_to_max(i%2,i)
def zero_to_max_alterate_step_instructions_non_shuffled():
    for i in range (0,N_neurons):
        _alternate_step_zero_to_max(0,i)


def _one_step(neuron:int,iteration:int, frequency:float):
    instruction_start_times[neuron,iteration]=iteration*interval_duration+neuron_offsets[neuron]
    instruction_end_times[neuron,iteration]=(iteration+1)*interval_duration+neuron_offsets[neuron]
    frequencies[neuron,iteration]=frequency
def _two_step(neuron:int,iteration:int, frequency:float):
    instruction_start_times[neuron,iteration]=iteration*interval_duration+neuron_offsets[neuron]
    instruction_end_times[neuron,iteration]=(iteration+1)*interval_duration+neuron_offsets[neuron]
    frequencies[neuron,iteration]=frequency
def _alternate_step_one_to_two(start:int,neuron, frequency1, frequency2):
    for j in range (0,total_time_intervals):
        if ((j+start)%2==0):
            _one_step(neuron,j, frequency1)
        else:
            _two_step(neuron,j, frequency2)
def one_to_two_alterate_step_instructions_shuffled(frequency1:float, frequency2:float):
    for i in range (0,N_neurons):
        _alternate_step_one_to_two(i%2,i, frequency1, frequency2)
def one_to_two_alterate_step_instructions_non_shuffled(frequency1:float, frequency2:float):
    for i in range (0,N_neurons):
        _alternate_step_one_to_two(0,i, frequency1, frequency2)

### Write instructions functions
def generate_instructions_from_arrays(neuron_ids:np.array, instruction_start_times:np.array, instruction_end_times:np.array,frequencies:np.array):
    with open (filename, 'w') as text_file:
        for i in range (0,N_neurons):
            for j in range (0,total_time_intervals):
                text_file.write('> '+ str(i)+' '+str(instruction_start_times[i,j])+' '+str(instruction_end_times[i,j])+ ' ' + str(frequencies[i,j]))
                if i+1 == N_neurons and j+1 == total_time_intervals:
                    pass
                else:
                    text_file.write('\n')
    text_file.close()


## Run the code
set_up_arrays_random_max_frequency(5.7, False)
one_to_two_alterate_step_instructions_non_shuffled(8.9, 4.6)
generate_instructions_from_arrays(neuron_ids, instruction_start_times, instruction_end_times, frequencies)