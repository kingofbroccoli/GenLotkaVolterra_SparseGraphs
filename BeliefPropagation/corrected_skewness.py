__author__ = 'david'

import numpy as np
from scipy.integrate import simpson
from scipy.special import erf
from scipy.optimize import fsolve
from scipy.stats import truncnorm
import matplotlib.pyplot as plt


def read_distribution(path, filename):
    fin = open(f'{path}/{filename}', 'r')
    nvals = []
    distr = []
    for line in fin:
        if line.startswith('#'):
            continue
        fields = line.split()
        if float(fields[0]) == 0:
            continue
        nvals.append(float(fields[0]))
        distr.append(float(fields[1]))
    return np.array(nvals), np.array(distr)


def compute_moments(nvals, distr):
    mean = simpson(nvals * distr, nvals)
    variance = simpson((nvals - mean) ** 2 * distr, nvals)
    skewness = simpson((nvals - mean) ** 3 * distr, nvals)
    kurtosis = simpson((nvals - mean) ** 4 * distr, nvals)
    return mean, variance, skewness, kurtosis


def alpha(mean_trunc, sigma_trunc, lower_bound):
    return (lower_bound - mean_trunc) / sigma_trunc


def residuals(vals, *params):
    mean, variance, lower_bound = params
    x0_trunc, sigma_trunc = vals
    a = alpha(x0_trunc, sigma_trunc, lower_bound)
    phi = np.exp(-0.5 * a ** 2) / np.sqrt(2 * np.pi)
    Phi = 0.5 * (1 + erf(a / np.sqrt(2)))
    mean_fit = x0_trunc + sigma_trunc * phi / (1 - Phi)
    variance_fit = sigma_trunc ** 2 * (1 + a * phi / (1 - Phi) - (phi / (1 - Phi)) ** 2)
    return [mean_fit - mean, variance_fit - variance]


def fit_truncated_gaussian(mean, variance, lower_bound):
    initial_guess = [mean, np.sqrt(variance)]
    params = (mean, variance, lower_bound)
    x0_trunc, sigma_trunc = fsolve(residuals, initial_guess, args=params)
    return x0_trunc, sigma_trunc


def print_diff_with_trunc_gaussian(filename, nvals, distr, x0_trunc, sigma_trunc, lower_bound):
    a = alpha(x0_trunc, sigma_trunc, lower_bound)
    trunc_gauss = truncnorm.pdf(nvals, a, np.inf, loc=x0_trunc, scale=sigma_trunc)
    rel_local_error = (distr - trunc_gauss) / trunc_gauss
    fout = open(filename, 'w')
    fout.write(f'# x0_trunc = {x0_trunc}\n')
    fout.write(f'# sigma_trunc = {sigma_trunc}\n')
    fout.write(f'# lower_bound = {lower_bound}\n')
    fout.write("# n  (P_BP-P_GAUSS)\n")
    for i in range(len(nvals)):
        fout.write(f'{nvals[i]} {rel_local_error[i]}\n')
    fout.close()


def print_trunc_gaussian(filename, x0_trunc, sigma_trunc, lower_bound, max_n, npoints=1000):
    a = alpha(x0_trunc, sigma_trunc, lower_bound)
    x = np.linspace(lower_bound, max_n, npoints)
    a = alpha(x0_trunc, sigma_trunc, lower_bound)
    trunc_gauss = truncnorm.pdf(x, a, np.inf, loc=x0_trunc, scale=sigma_trunc)
    fout = open(filename, 'w')
    fout.write(f'# x0_trunc = {x0_trunc}\n')
    fout.write(f'# sigma_trunc = {sigma_trunc}\n')
    fout.write(f'# lower_bound = {lower_bound}\n')
    fout.write("# n  P_GAUSS(n)\n")
    for i in range(len(x)):
        fout.write(f'{x[i]} {trunc_gauss[i]}\n')
    fout.close()
    


def plot_distribution(nvals, distr, x0_trunc, sigma_trunc, lower_bound, mu, T):
    plt.plot(nvals, distr, label=f'T={T} mu={mu}', marker='x', linestyle='None')
    x = np.linspace(lower_bound, nvals[-1], 1000)
    a = alpha(x0_trunc, sigma_trunc, lower_bound)
    trunc_gauss = truncnorm.pdf(x, a, np.inf, loc=x0_trunc, scale=sigma_trunc)
    plt.plot(x, trunc_gauss, linestyle='--', color='black')
    plt.xlabel('n')
    plt.ylabel('P(n)')
    plt.yscale('log')
    plt.ylim(np.min(distr[distr > 0]), np.max(distr) * 2)
    plt.legend()
    plt.show()


def plot_difference(nvals, distr, x0_trunc, sigma_trunc, lower_bound, mu, T):
    a = alpha(x0_trunc, sigma_trunc, lower_bound)
    trunc_gauss = truncnorm.pdf(nvals, a, np.inf, loc=x0_trunc, scale=sigma_trunc)
    rel_local_error = (distr - trunc_gauss) / trunc_gauss
    plt.plot(nvals, rel_local_error, label=f'T={T} mu={mu}', linestyle='-', color='black')
    plt.xlabel('n')
    plt.ylabel(r'$P_{BP}(n)-P_{GAUSS}(n)$')
    plt.xlim(0.1, 1.5)
    plt.legend()
    plt.show()



def main():
    N = "128"
    # T = "0.050"
    # mu_list = ["0.250", "0.260", "0.270", "0.280", "0.290", "0.300", "0.310", "0.320", "0.330"]
    # T = "0.036"
    # mu_list = ["0.030", "0.040", "0.050", "0.060", "0.070", "0.080", "0.090"]
    T = "0.030"
    # mu_list = ["0.020", "0.040", "0.050", "0.060", "0.070", "0.080", "0.090", "0.100", "0.110", "0.120", "0.130"]
    mu_list = ["0.020"]
    # T = "0.025"
    # mu_list = ["0.060", "0.070", "0.080", "0.090", "0.100", "0.110", "0.120", "0.130"]
    # mu_list = ["0.060", "0.130"]
    # T = "0.010"
    # mu_list = ["0.190", "0.200", "0.210", "0.220", "0.230"]
    # T = "0.001"
    # mu_list = ["0.100", "0.200", "0.300", "0.310", "0.320", "0.330", "0.340"]
    
    T_float = float(T)
    ndigits = 6

    path = '/media/david/Data/UH/Grupo_de_investigacion/Ecology/BP/Code/Graph/Debug'
    # path = '/media/david/Seagate Expansion Drive/Salva/Salva_Data_Investigacion/Grupo_de_investigacion/Ecology/BP/Code/Graph/Debug'
    
    for mu in mu_list:
        filename = f'BP_cont_LV_sigma0_true_marginals_gaussian_part_N_{N}_T_{T}_mu_{mu}.txt'
        nvals, distr = read_distribution(path, filename)
        lower_bound = nvals[0]
        mean, variance, skewness, kurtosis = compute_moments(nvals, distr)
        x0_trunc, sigma_trunc = fit_truncated_gaussian(mean, variance, lower_bound)
        a = alpha(x0_trunc, sigma_trunc, lower_bound)
        mean_trunc, variance_trunc, skewness_trunc, kurtosis_trunc = \
            truncnorm.stats(a, np.inf, loc=x0_trunc, scale=sigma_trunc, moments='mvsk')
        
        print(f'{T}  {mu}  {skewness:.{ndigits}f}  {skewness_trunc:.{ndigits}f}  {kurtosis:.{ndigits}f}  {kurtosis_trunc:.{ndigits}f}  {mean_trunc:.{ndigits}f}  {(variance_trunc / T_float):.{ndigits}f}  {x0_trunc:.{ndigits}f}  {(sigma_trunc**2 / T_float):.{ndigits}f}')
        filename_diff = f'BP_cont_LV_sigma0_diff_true_marginals_fit_trunc_gauss_N_{N}_T_{T}_mu_{mu}.txt'
        print_diff_with_trunc_gaussian(filename_diff, nvals, distr, x0_trunc, sigma_trunc, lower_bound)
        filename_gauss = f'BP_cont_LV_sigma0_true_marginals_fit_trunc_gauss_N_{N}_T_{T}_mu_{mu}.txt'
        print_trunc_gaussian(filename_gauss, x0_trunc, sigma_trunc, lower_bound, nvals[-1])
        # plot_distribution(nvals, distr, x0_trunc, sigma_trunc, lower_bound, mu, T)
        # plot_difference(nvals, distr, x0_trunc, sigma_trunc, lower_bound, mu, T)
    

    return 0


if __name__ == '__main__':
    main()
