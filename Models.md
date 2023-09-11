# Models
## MongilloSynapse
Pierre
## MongilloSynapseContinuous
Pierre
## PRGSynapse
Pierre
## PowerLawSynapse
Pierre
## AdjacencyMatrixConnectivity
Ulzii
(Now it can contain multiple synapses)
## PoissonConnectivity
This connectivity class is theoretically equivalent to RandomConnectivity, but instead of getting a set of synapses from the probability of connection and the source neurons which is fixed, you sample a poisson distribution with the class' parameter `connectivity_connectionLambda` as the distribution's lambda. This allows for multiple synapses between two neurons with a random process.
## TraceRBranchedHSTDP
*(Trace based heterosynaptic plasticity model with multiple resource-limited dendritic branches, by Antoni Bertolin Monferrer)*

This model is a discrete implementation of the conjunction of the heterosynaptic expansion of trace-based STDP with a resource limitation per branch. The model can be used on any type of neuron with the parameter `synapses_N_N_pmodel_type TraceRBranchedHSTDP`, and the basic morphological parameters are `synapses_N_N_pmodel_branchLength`, `synapses_N_N_pmodel_synapticGap`, `synapses_N_N_pmodel_dendriteBranchings`, and `synapses_N_N_pmodel_synapseAllocation`. The synaptic gap and branch length define the number of available synaptic slots and the distance between them. This two are in units of micrometres, and all spatial constants in the models are expressed in the same unit. The branchings parameter does not define the total branching points in the dendritic tree, but the number of branching points that each path in the dendritic tree goes through. The total amount of branching points in the tree is given by $2^{branchings}$ . The synapse allocation refers to the order in which synapses are allocated to the synaptic branch. Once the branch has been allocated (this pertains to the Synapse object), the synapses can be allocated in that branch in an ordered manner (from position 0 to last position) or randomly (a random unoccupied position is chosen).

The model's behaviour can be separated into two parts: the traces and the definition of synaptic weight through resource limitation.

### Traces

In this model, every synaptic spine contains three traces (presynaptic trace, cooperativity trace, alpha trace) plus its weight, and the neuron as a whole has a trace for the postsynaptic spikes. 
#### STDP traces
Each spine's presynaptic trace and the postsynaptic trace are defined by the same equation:
$$\frac{dT_{i/p}}{dt}=-\frac{T_{i/p}}{\tau_{t}}+S_{i/p}(t)$$

$$S_{i}(t)= \sum_{spike}\delta(t-t_{i,spike})$$
     
Where $\tau_{t}$ is the decay of the trace in seconds, defined by `synapses_N_N_pmodel_tauSTDP`.
Both of which aim to reproduce the behaviour of the eligibility trace implementation of STDP through their interaction, while taking into account the difference in effect between a presynaptic/postsynaptic spike train vs a single spike in the weight. 
#### Cooperativity traces
Each spine's cooperativity trace is defined as:
$$\frac{dC_{i}}{dt}=-\frac{C_{i}}{\tau_{c}}+\sum_{j=0, i\neq j}^{N}(S_{j}(t)exp^{(\frac{-d_{i,j}}{\lambda})})$$
Where $S_{j}(t)$ is either zero or one depending on whether there was a presynaptic spike in the presynaptic trace of synapse $j$, $d_{i,j}$ is the distance between synapses $i$ and $j$. $\tau_{t}$ is the decay of the trace in seconds, defined by `synapses_N_N_pmodel_coopTau`, and $\lambda$ is the decay constant that defines the spatial profile of synaptic cooperativity in the model, defined by `synapses_N_N_pmodel_coopProfile`.
Which aim to reproduce heterosynaptic plasticity behaviour. There is an assumption of instantaneous effect on the cooperativity trace. This assumption is made for practical reasons.
#### Alpha (resource) traces
Each spine's alpha (or resource) trace is defined as:

$$ \frac{d\alpha_{i}}{dt}=-\frac{\alpha_{i}-\alpha_{basal}}{\tau_{\alpha}}+F_i(t)$$



$$F_i(t)= \begin{cases} \alpha_{st}T_{i}(1+C_i) &\quad S_p(t)=1 \newline -\alpha_{st}T_{i}(1-C_i) &\quad S_p(t)=0 \wedge S_i(t)=1 \newline 0 &\quad S_p(t)=0 \wedge S_i(t)=0 \end{cases}$$


where $\alpha_{basal}$ is the basal value, defined by `synapses_N_N_pmodel_basalAlpha`. $\alpha_{st}$ is the constant amount by which the trace increases or decreases per event, defined by `synapses_N_N_pmodel_baseAlphaIncrease`. the decay constant $\tau_{\alpha}$ is the trace decay to the basal value, defined by `synapses_N_N_pmodel_alphaTau`.
These traces represent the spine's access to resources in the neuron, and they will be used in the resource normalization of weights. Changes in these traces are not capped upwards but are not allowed to become negative, as the excitatory-inhibitory inversion is undesirable in this model

### Resource normalization of weights
Weights in each dendritic branch are normalized in this model according to the total amount of trace in the branch plus a constant. The definition of weight is:

$$w_{i}(t)=\beta\frac{\alpha_{trace, i}(t)}{\omega + \sum_i (\alpha_{trace, i}(t))}$$

where $\beta$ is the total amount of weight in the branch (in dmV/spike) to be distributed according to the fraction and is defined by `synapses_N_N_pmodel_betaResourcePool`, $\omega$ represents the amount of resources that are constantly contained in the branch outside the spines (adimensional) and is defined by `synapses_N_N_pmodel_omegaOffset`. $w_{i}(t)$ is the weight of the i-th spine at time $t$$t$, and the sum of the N $\alpha_{trace}(t)$ is the total amount of traces at time $t$. 
This resource limiting normalization aims to reproduce another heterosynaptic effect, that is the effect of allocating available resources unequally along the branch. This normalization also fulfills the need for capping the synaptic weights in a more organic way, as the total weight of a dendritic branch will never go beyond $\beta$.

## MonoDendriteSTDP models
The models within this tree of classes work with the $\theta$ trace, also called cooperativity. The models were designed by Saif Ahmed and were integrated with the current code by Antoni Bertolin. These models also use dendritic length and synaptic gap for synapse allocation (see Test 10 in [Tests](README_Tests.md)).All models calculate every timestep said trace, which is determined by the equation:

$$ \frac{d\theta_i}{dt}=-\frac{\theta_i}{\tau_{\theta}}+\sum_{\forall k\neq i}w_iw_k\sum_{\forall j}\delta(t-t_{k,pre}^{(j)})\exp(-\frac{|x_i-x_k|}{\lambda_{dist}})\exp(-\frac{|t_{i,pre}^{(last)}-t_{k,pre}^{(j)}|}{\tau_{delay}})+\\
$$

$$
 \quad \sum_{\forall k\neq i}w_iw_k\sum_{\forall j}\delta(t-t_{i,pre}^{(j)})\exp(-\frac{|x_i-x_k|}{\lambda_{dist}})\exp(-\frac{|t_{i,pre}^{(j)}-t_{k,pre}^{(last)}|}{\tau_{delay}})
$$

where $\theta_i$ represents synaptic cooperativity, $\tau_{\theta}$(`synapses_N_N_pmodel_heterosynapticThetaDecay`) its decay, $t_{k,pre}^{(j)}$ the arrival time of the j-th presynaptic input at synapse k (when it is $last$ it refers to the closest input in time in the past). $x_i$ represents the distance of synapse i relative to the neuronal cell body or soma. $\lambda_{dist}$(`synapses_N_N_pmodel_intersynapseDistanceDecay`) is the decay constant of said distance. $\tau_{delay}$(`synapses_N_N_pmodel_intersynapseSpikeTimingDecay`) represents the decay constant of the temporal decay that weakens the cooperative interactions.

Then the weight is determined by the equations:

$$\frac{dw_i}{dt}=H_{LTP}(\theta_i(t))A_{LTP}\sum_{\forall j}(\delta(t-t_{post}^{(j)})K(t_{post}^{(j)}-t_{i,pre}^{(last)})) + \ 
\newline H_{LTD}(\theta_i(t))A_{LTD}\sum_{\forall j}(\delta(t-t_{i,pre}^{(last)})K(t_{i,pre}^{(j)}-t_{post}^{(last)})) 
 $$
 
 $$
H_{LTP}(\theta)= B_{LTP}+I_{LTP}(1-\exp(-\alpha\theta))
 $$
 
$$
H_{LTD}(\theta)= B_{LTD}-D_{LTD}(1-\exp(-\beta\theta))
$$

Where $A_{LTP/LTD}$ (`synapses_N_N_pmodel_preFactorLtp, synapses_N_N_pmodel_preFactorLtd`) represents the maximal potentiation/depression that a synapse can undergo in one pair-based event. $K$ (depends on the class) represents the kernel function used for potentiation and depression, and $w_i$ represents the weight of the $i$-th spine. $H_{LTP/LTD}$ represents the heterosynaptic effects on LTP and LTD, defined above. $B_{LTP/LTD}$ (`synapses_N_N_pmodel_baseLtp,synapses_N_N_pmodel_baseLtd`) represents the baseline LTP/LTD with no heterosynaptic cooperativity. $I_{LTP}$ (`synapses_N_N_pmodel_incrLtp`) is the base increase in LTP due to cooperativity, which is modulated by $\alpha$ (`synapses_N_N_pmodel_alphaLtp`). $D_{LTD}$ (`synapses_N_N_pmodel_decrLtd`) is the base decrease of LTD due to cooperativity, modulated by $\beta$ (`synapses_N_N_pmodel_betaLtd`). For time notation, look at the previous paragraph.

Every class inside the MonoDendriteSTDP family follows these equations, and each model differs in the kernel $K$, and has particular parameters that define said kernel.
