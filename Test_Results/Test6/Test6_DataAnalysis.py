#Translation from MATLAB data analysis to python. Currently not finished
import numpy as np
from scipy.optimize import fsolve, least_squares
from scipy.special import erf, zeta
from scipy.integrate import quad
import matplotlib.pyplot as plt
import os

def Phi(mu, sigma, VR, tau_m, tau_s=0):
    def Phi_Integrand(z):
        mark = (z < -4)
        y = np.zeros_like(z)
        y[~mark] = Phi_Integrand_Full(z[~mark])
        y[mark] = Phi_Integrand_Expansion(z[mark])
        return y
    
    def Phi_Integrand_Full(z):
        return np.exp(z ** 2) * (1 + erf(z))
    
    def Phi_Integrand_Expansion(z):
        if z > 0:
            raise ValueError('Error')
        
        z = -z
        zt = 2 * z ** 2
        return 1 / np.sqrt(np.pi) / z * (1 - 1 / zt + 3 / zt ** 2 - 3 * 5 / zt ** 3)
    if len(sigma) == 1:
        sigma = np.ones_like(mu) * sigma
    
    alpha = np.sqrt(2) * np.abs(zeta(0.5))
    d_tau_s = np.sqrt(tau_s / tau_m) * alpha / 2

    z = np.zeros_like(mu)
    for i in range(len(mu)):
        y1 = (VR - tau_m * mu[i]) / (sigma[i] * np.sqrt(tau_m)) + d_tau_s
        y2 = (1 - tau_m * mu[i]) / (sigma[i] * np.sqrt(tau_m)) + d_tau_s

        if Phi_Integrand(y2) > 1e10:
            z[i] = 0
        else:
            z[i] = (tau_m * np.sqrt(np.pi) * quad(Phi_Integrand, y1, y2)) ** (-1)

    return z

def import_parameters(filename):
    params = {}

    with open(filename, 'r') as file:
        for line in file:
            arr = line.strip().split()
            if len(arr) == 0:
                continue
            if arr[0] == 'neurons_0_noNeurons':
                N_E = int(arr[1])
                params['N_E'] = N_E
            elif arr[0] == 'neurons_1_noNeurons':
                N_I = int(arr[1])
                params['N_I'] = N_I
            elif arr[0] == 'stimulus_meanCurrent':
                params['mu_E_ext'] = float(arr[1])
                params['mu_I_ext'] = float(arr[2])
            elif arr[0] == 'stimulus_sigmaCurrent':
                params['mu_E_ext_sigma'] = float(arr[1])
                params['mu_I_ext_sigma'] = float(arr[2])
            elif arr[0] == 'dt':
                params['dt'] = float(arr[1])
            elif arr[0] == 'synapses_0_0_J':
                params['J_EE'] = float(arr[1])
            elif arr[0] == 'synapses_0_1_J':
                params['J_EI'] = float(arr[1])
            elif arr[0] == 'synapses_1_0_J':
                params['J_IE'] = float(arr[1])
            elif arr[0] == 'synapses_1_1_J':
                params['J_II'] = float(arr[1])
            elif arr[0] == 'neurons_0_vReset':
                params['V_R'] = float(arr[1])
            elif arr[0] == 'neurons_0_tauM':
                params['tau_m'] = float(arr[1])
            elif arr[0] == 'synapses_0_0_connectivity_ConnectionProba':
                params['conn_0_0'] = float(arr[1])

    params['N'] = params['N_E'] + params['N_I']

    return params

def import_data_file(filenames, start_row=26, end_row=-1):
    delimiter = '\t'

    # Format string for each line of text
    #format_spec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]'
    # Read the text file
    filename=[filename for filename in filenames if filename.endswith("Data.dat")]
    with open(filename[0], 'r') as file:
        file_data = file.readlines()

    # Extract the relevant data rows
    data_rows = file_data[start_row - 1:end_row]

    # Initialize empty arrays to store the data
    t = []
    v_e = []
    v_i = []
    nu_e = []
    nu_i = []
    mu_ext_e = []
    mu_ext_i = []
    mu_e = []
    mu_i = []
    sigma_e = []
    sigma_i = []

    # Process each data row
    for row in data_rows:
        line_data = row.strip().split(delimiter)
        t.append(float(line_data[0]))
        v_e.append(float(line_data[1]))
        v_i.append(float(line_data[2]))
        nu_e.append(float(line_data[3]))
        nu_i.append(float(line_data[4]))
        mu_ext_e.append(float(line_data[5]))
        mu_ext_i.append(float(line_data[6]))
        mu_e.append(float(line_data[11]))
        mu_i.append(float(line_data[12]))
        sigma_e.append(float(line_data[15]))
        sigma_i.append(float(line_data[16]))

    # Create a dictionary to store the data
    data = {
        't': t,
        'V': [v_e, v_i],
        'nu': [nu_e, nu_i],
        'mu_Ext': [mu_ext_e, mu_ext_i],
        'mu': [mu_e, mu_i],
        'sigma': [sigma_e, sigma_i]
    }

    # Import the parameters
    filename_parameters=[filename for filename in filenames if filename.endswith("Parameters.txt")]
    data['params'] = import_parameters(filename_parameters[0])

    return data

def GetInvPhi(nu, sigma, VR, tau_m, mu_guess=None):
    opts = {'disp': False, 'xtol': 1e-8}

    if mu_guess is None:
        mu_guess = 1 / tau_m

    mu = np.zeros_like(nu)

    for ii in range(len(nu)):
        nu_ = nu[ii]
        sigma_ = sigma[ii]

        if ii > 0:
            mu[ii] = fsolve(lambda x: f(x, nu_, sigma_, VR, tau_m), mu[ii-1], **opts)
        else:
            mu[ii] = fsolve(lambda x: f(x, nu_, sigma_, VR, tau_m), mu_guess, **opts)
    def f(mu, nu, sigma, VR, tau_m):
        return Phi(mu, sigma, VR, tau_m) - nu
    return mu

def ComputeSteadyStates(params):
    # Numerical parameters
    J_EE = abs(params['J_EE'])
    J_EI = abs(params['J_EI'])
    J_IE = abs(params['J_IE'])
    J_II = abs(params['J_II'])
    N = params['N']
    VR = params['V_R']
    tau_m = params['tau_m']
    mu_E_ext = params['mu_E_ext']
    mu_I_ext = params['mu_I_ext']
    mu_E_ext_sigma = params['mu_E_ext_sigma']
    mu_I_ext_sigma = params['mu_I_ext_sigma']
    N_E = params['N_E']
    N_I = params['N_I']
    conn_prob = params['conn_0_0']
    C_E = conn_prob * N_E  # Number of excitatory inputs per neuron
    C_I = conn_prob * N_I  # Number of inhibitory inputs per neuron
    opts = 'N_scaling'  # 'C_scaling'
    
    # Rescale weights
    if opts == 'N_scaling':
        J_EE_rescaled = J_EE * C_E / N
        J_EI_rescaled = J_EI * C_I / N
        J_II_rescaled = J_II * C_I / N
        J_IE_rescaled = J_IE * C_E / N
        NC = N
    elif opts == 'C_scaling':
        J_EE_rescaled = J_EE
        J_EI_rescaled = J_EI * np.sqrt(C_I / C_E)
        J_II_rescaled = J_II * np.sqrt(C_I / C_E)
        J_IE_rescaled = J_IE
        NC = C_E
        
    J = np.array([[J_EE_rescaled, -J_EI_rescaled], [J_IE_rescaled, -J_II_rescaled]])

    z = fsolve(balanced_state, [1, 1])
    nu_E_balanced = z[0]
    nu_I_balanced = z[1]
    
    sigma_E_balanced, sigma_I_balanced = get_sigma(nu_E_balanced, nu_I_balanced)
    mu_E_balanced = GetInvPhi(nu_E_balanced, sigma_E_balanced, VR, tau_m)
    mu_I_balanced = GetInvPhi(nu_I_balanced, sigma_I_balanced, VR, tau_m)
    
    x = np.array([nu_E_balanced, nu_I_balanced, sigma_E_balanced, sigma_I_balanced])
    z = least_squares(BalancedStateCondition_Finite_C, x, bounds=([0, 0, 0, 0], [1000, 1000, 1000, 1000]))
    nu_E_finite = z[0]
    nu_I_finite = z[1]
    sigma_E_finite = z[2]
    sigma_I_finite = z[3]

    mu_E_finite = GetInvPhi(nu_E_finite, sigma_E_finite, VR, tau_m, mu_E_balanced)
    mu_I_finite = GetInvPhi(nu_I_finite, sigma_I_finite, VR, tau_m, mu_I_balanced)

    print('************************************************')
    print('Condition 1      : mu_I_ext/mu_E_ext < J_II/J_EI < J_IE/J_EE')
    print('for this network :', mu_I_ext / mu_E_ext, '<', J_II_rescaled / J_EI_rescaled, '<', J_IE_rescaled / J_EE_rescaled)
    if (mu_I_ext / mu_E_ext) > (J_II_rescaled / J_EI_rescaled) or (mu_I_ext / mu_E_ext) > (J_IE_rescaled / J_EE_rescaled) or (J_II_rescaled / J_EI_rescaled) > (J_IE_rescaled / J_EE_rescaled):
        print('Warning: conditions not met!')
    print('Condition 2     : J_IE > J_EE')
    print('for this network :', J_IE_rescaled, '>', J_EE_rescaled)
    if J_EI_rescaled < J_EE_rescaled:
        print('Warning: J_EI < J_EE')
    print('det(J)  =', np.linalg.det(J))
    print('************************************************')

    print('Value   = finite C -- balanced state')
    print('nu_E    =', nu_E_finite, '  -- ', nu_E_balanced)
    print('nu_I    =', nu_I_finite, '  -- ', nu_I_balanced)
    print('sigma_E =', sigma_E_finite, '  -- ', sigma_E_balanced)
    print('sigma_I =', sigma_I_finite, '  -- ', sigma_I_balanced)
    print('mu_E    =', mu_E_finite, '  -- ', mu_E_balanced)
    print('mu_I    =', mu_I_finite, '  -- ', mu_I_balanced)

    def balanced_state(x):
        nuE_ = x[0]
        nuI_ = x[1]
        z = np.zeros(2)
        z[0] = mu_E_ext + J_EE_rescaled * nuE_ - J_EI_rescaled * nuI_
        z[1] = mu_I_ext + J_IE_rescaled * nuE_ - J_II_rescaled * nuI_
        return z
    def get_sigma(nu_E, nu_I):
        if opts == 'N_scaling':
            sE = np.sqrt(N / C_E * J_EE_rescaled ** 2 * nu_E + N / C_I * J_EI_rescaled ** 2 * nu_I + mu_E_ext_sigma ** 2)
            sI = np.sqrt(N / C_E * J_IE_rescaled ** 2 * nu_E + N / C_I * J_II_rescaled ** 2 * nu_I + mu_I_ext_sigma ** 2)
        elif opts == 'C_scaling':
            sE = np.sqrt(J_EE_rescaled^2*nu_E + C_E/C_I*J_EI_rescaled^2*nu_I + mu_E_ext_sigma^2)
            sI = np.sqrt(J_IE_rescaled^2*nu_E + C_E/C_I*J_II_rescaled^2*nu_I + mu_I_ext_sigma^2)
        return sE, sI
    def BalancedStateCondition_Finite_C(x):
        nu_E = x[0]
        nu_I = x[1]
        sigma_E = x[2]
        sigma_I = x[3]

        mu_E = GetInvPhi(nu_E, sigma_E, VR, tau_m, mu_E_balanced)  # Eq. (14-E) Suppl.
        mu_I = GetInvPhi(nu_I, sigma_I, VR, tau_m, mu_I_balanced)  # Eq. (14-I) Suppl.

        z = np.zeros(4)
        z[0] = mu_E_ext + J_EE_rescaled * nu_E - J_EI_rescaled * nu_I - mu_E / np.sqrt(NC)
        z[1] = mu_I_ext + J_IE_rescaled * nu_E - J_II_rescaled * nu_I - mu_I / np.sqrt(NC)  # Eq. (25) Suppl.
        sE, sI = get_sigma(nu_E, nu_I)
        z[2] = sigma_E - sE  # Eq. (11-E) Suppl.
        z[3] = sigma_I - sI  # Eq. (11-I) Suppl.

        return z
    return nu_E_balanced, nu_I_balanced, mu_E_balanced, mu_I_balanced, sigma_E_balanced, sigma_I_balanced

def ComputeBalancedState(params):
    J_EE = abs(params.J_EE)
    J_EI = abs(params.J_EI)
    J_IE = abs(params.J_IE)
    J_II = abs(params.J_II)

    N = params.N
    VR = params.V_R
    tau_m = params.tau_m

    mu_E_ext = params.mu_E_ext
    mu_I_ext = params.mu_I_ext

    mu_E_ext_sigma = params.mu_E_ext_sigma
    mu_I_ext_sigma = params.mu_I_ext_sigma

    N_E = params.N_E
    N_I = params.N_I

    connProb = params.conn_0_0

    C_E = connProb * N_E
    C_I = connProb * N_I

    opts = 'N_scaling'

    if opts == 'N_scaling':
        J_EE_rescaled = J_EE * C_E / N
        J_EI_rescaled = J_EI * C_I / N
        J_II_rescaled = J_II * C_I / N
        J_IE_rescaled = J_IE * C_E / N
    elif opts == 'C_scaling':
        J_EE_rescaled = J_EE
        J_EI_rescaled = J_EI * np.sqrt(C_I / C_E)
        J_II_rescaled = J_II * np.sqrt(C_I / C_E)
        J_IE_rescaled = J_IE

    J = np.array([[J_EE_rescaled, -J_EI_rescaled], [J_IE_rescaled, -J_II_rescaled]])

    def BalancedState(x):
        nuE_ = x[0]
        nuI_ = x[1]

        return [
            mu_E_ext + J_EE_rescaled * nuE_ - J_EI_rescaled * nuI_,
            mu_I_ext + J_IE_rescaled * nuE_ - J_II_rescaled * nuI_
        ]

    z = fsolve(BalancedState, [1, 1])
    nu_E_balanced = z[0]
    nu_I_balanced = z[1]

    sigma_E_balanced, sigma_I_balanced = GetSigma(nu_E_balanced, nu_I_balanced)
    mu_E_balanced = GetInvPhi(nu_E_balanced, sigma_E_balanced, VR, tau_m)
    mu_I_balanced = GetInvPhi(nu_I_balanced, sigma_I_balanced, VR, tau_m)

    def GetSigma(nu_E, nu_I):
        if opts == 'N_scaling':
            sE = np.sqrt(N / C_E * J_EE_rescaled ** 2 * nu_E + N / C_I * J_EI_rescaled ** 2 * nu_I + mu_E_ext_sigma ** 2)
            sI = np.sqrt(N / C_E * J_IE_rescaled ** 2 * nu_E + N / C_I * J_II_rescaled ** 2 * nu_I + mu_I_ext_sigma ** 2)
        elif opts == 'C_scaling':
            sE = np.sqrt(J_EE_rescaled ** 2 * nu_E + C_E / C_I * J_EI_rescaled ** 2 * nu_I + mu_E_ext_sigma ** 2)
            sI = np.sqrt(J_IE_rescaled ** 2 * nu_E + C_E / C_I * J_II_rescaled ** 2 * nu_I + mu_I_ext_sigma ** 2)
        return sE, sI

    return nu_E_balanced, nu_I_balanced, mu_E_balanced, mu_I_balanced, sigma_E_balanced
data = []
fullData = {}

# ***********************************
# Import data
# ***********************************
# TODO: Set here folder in which Test6_.. folders are located
base_folder = os.getcwd()
# ***********************************

folders = folders = [f.name for f in os.scandir(base_folder) if f.is_dir()]
files=[]
for folder in folders:
    if folder.startswith('Test6') and os.path.isdir(folder):
        for file in os.listdir(folder):
            if file.endswith("Data.dat") or file.endswith("Parameters.txt"):
                dat_name = os.path.join(base_folder, folder, file)
                files.append(dat_name)
        print(f'Loading {folder}..')
        data.append(import_data_file(files))
# Clear previous variables
# fullData = {}

# # Import data
# base_folder = '../'

# folders = [f.name for f in os.scandir(base_folder) if f.is_dir()]
# data = []

# for folder_name in folders:
#     if len(folder_name) >= 4 and folder_name[:4] == 'Test':
#         dat_name = next(os.scandir(os.path.join(base_folder, folder_name)), None)
#         if dat_name is not None and dat_name.name.endswith('.dat'):
#             fullname = dat_name.path
#             print(f'Loading {fullname}..')
#             data.append(import_data_file(fullname))  # Import_DataFile needs to be defined

# Data processing
fullData = {key: [] for key in ['N', 'nu_E', 'nu_I', 'mu_E', 'mu_I', 'sigma_E', 'sigma_I', 'nu_E_Eq3', 'nu_I_Eq3', 'nu_E_finite', 'nu_I_finite', 'mu_E_finite', 'mu_I_finite', 'sigma_E_finite', 'sigma_I_finite']}

for datum in data:
    fullData['N'].append(datum['params']['N'])
    fullData['nu_E'].append(datum['nu'][0][-1])
    fullData['nu_I'].append(datum['nu'][1][-1])

    fullData['mu_E'].append(datum['mu'][0][-1])
    fullData['mu_I'].append(datum['mu'][1][-1])

    fullData['sigma_E'].append(datum['sigma'][0][-1])
    fullData['sigma_I'].append(datum['sigma'][1][-1])

    fullData['nu_E_Eq3'].append(Phi(np.full(1,datum['mu'][0][-1]), np.full(1,datum['sigma'][0][-1]), 0.0, 0.01))  # Phi needs to be defined
    fullData['nu_I_Eq3'].append(Phi(np.full(1,datum['mu'][1][-1]), np.full(1,datum['sigma'][1][-1]), 0.0, 0.01))  # Phi needs to be defined

    steady_states = ComputeSteadyStates(datum['params'])  # ComputeSteadyStates needs to be defined

    fullData['nu_E_finite'].append(steady_states[0])
    fullData['nu_I_finite'].append(steady_states[1])
    fullData['mu_E_finite'].append(steady_states[2])
    fullData['mu_I_finite'].append(steady_states[3])
    fullData['sigma_E_finite'].append(steady_states[4])
    fullData['sigma_I_finite'].append(steady_states[5])

balanced_state = ComputeBalancedState(data[-1]['params'])  # ComputeBalancedState needs to be defined
O = np.ones(len(fullData['N']))

# Plotting
fig, axs = plt.subplots(1, 3)

color_1 = 'r'
color_2 = 'b'

axs[0].plot(fullData['N'], fullData['nu_E'], 'o', color=color_1, label='exc')
axs[0].plot(fullData['N'], fullData['nu_I'], 'o', color=color_2, label='inh')

axs[0].plot(fullData['N'], fullData['nu_E_finite'], '--', color=color_1)
axs[0].plot(fullData['N'], fullData['nu_I_finite'], '--', color=color_2)

axs[0].plot(fullData['N'], balanced_state[0]*O, '-.', color=color_1)
axs[0].plot(fullData['N'], balanced_state[1]*O, '-.', color=color_2)

axs[0].set_xlabel('N')
axs[0].set_ylabel('Firing rate (Hz)')

axs[1].plot(fullData['N'], fullData['mu_E'], 'o', color=color_1)
axs[1].plot(fullData['N'], fullData['mu_I'], 'o', color=color_2)

axs[1].plot(fullData['N'], fullData['mu_E_finite'], '--', color=color_1)
axs[1].plot(fullData['N'], fullData['mu_I_finite'], '--', color=color_2)

axs[1].plot(fullData['N'], balanced_state[2]*O, '-.', color=color_1)
axs[1].plot(fullData['N'], balanced_state[3]*O, '-.', color=color_2)

axs[1].set_xlabel('N')
axs[1].set_ylabel('Current $\\mu$')

axs[2].plot(fullData['N'], fullData['sigma_E'], 'o', color=color_1)
axs[2].plot(fullData['N'], fullData['sigma_I'], 'o', color=color_2)

axs[2].plot(fullData['N'], fullData['sigma_E_finite'], '--', color=color_1)
axs[2].plot(fullData['N'], fullData['sigma_I_finite'], '--', color=color_2)

axs[2].plot(fullData['N'], balanced_state[4]*O, '-.', color=color_1)
axs[2].plot(fullData['N'], balanced_state[5]*O, '-.', color=color_2)

axs[2].set_xlabel('N')
axs[2].set_ylabel('$\\sigma$')

fig.legend(loc='upper right', bbox_to_anchor=(1.2, 1), borderaxespad=0., labels=['exc', 'inh', 'simulation', 'finite', 'balanced'])
fig.tight_layout()
plt.savefig('N_convergence.pdf', format='pdf', bbox_inches='tight')
plt.show()
