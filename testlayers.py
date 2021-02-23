import numpy as np
import matplotlib.pyplot as plt
import tmm

degree = np.pi/180

# list of layer thicknesses in nm
d_list = [tmm.inf, 5, tmm.inf]
# list of refractive indices
n_list = [1, 2.68, 1]
nabs_list = [1 ,2.68+0.01j, 1]
# list of wavenumbers to plot in nm^-1
lams = np.linspace(1, 1.25, num=400)
# initialize lists of y-values to plot
Rnorm = []
Rabs = []
for lam in lams:
    # For normal incidence, s and p polarizations are identical.
    # I arbitrarily decided to use 's'.
    Rnorm.append(tmm.coh_tmm('s', n_list, d_list, 0, lam)['R'])
    #R45.append(tmm.unpolarized_RT(n_list, d_list, 45*degree, lam)['R'])
    Rabs.append(tmm.coh_tmm('s', nabs_list, d_list, 0, lam)['R'])

plt.figure()
plt.plot(lams, Rnorm, 'blue', lams, Rabs, 'purple')
plt.xlabel('$\lambda$ ($\mu$m)')
plt.ylabel('Fraction reflected')
plt.title('Modest page 56')
plt.show()

