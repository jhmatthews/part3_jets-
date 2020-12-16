from scipy.optimize import fsolve
import numpy as np 
import matplotlib.pyplot as plt 

# median jet power in erg/s
Qmedian = 1e45


# read data from file 
sigma = 1.5
seed = 38
time, flux = np.genfromtxt("Lightcurve_sig{:.1f}_seed{:d}.dat".format(sigma, seed), unpack=True)

# constants and so on 
C = 2.997925e10
unit_density = 6e-27  # this is the unit density in the code
PARSEC = 3.086E18
rho_j = 1e-4
width = 1.0 * PARSEC * 1000.0 # 1kpc jet nozzle width
area = np.pi  * (width)**2

# array to store the gammas 
gammas = np.zeros_like(flux)


for i in range(len(flux)):
    f = flux[i] * Qmedian

    # this is power over area over rho c**3
    X = f / area / rho_j / unit_density / (C ** 3) 
    func = lambda gmm: np.sqrt(1-(gmm**-2)) * ((gmm**2) - gmm) - X

    # find gamma 
    gammas[i] = fsolve(func, 3.0)[0]
    

data_to_save = np.column_stack((time, gammas))
fname = "gmm_Q{:.0f}_sig{:.1f}_eta{:.0f}_seed{:d}.dat".format(np.log10(Qmedian), sigma, -np.log10(rho_j), seed)
#fname = "data/Lightcurve_sig{:.1f}_seed{:d}.dat".format(sigma, seed)
print (fname)
np.savetxt(fname, data_to_save)
# make a plot
# plt.plot(time, gammas)
# plt.xlabel("t (Myr)")
# plt.ylabel(r"$\Gamma$")  
# plt.show()