import photovoltaic.semiconductors as semi
import photovoltaic.common as pv

print('Thermal voltage 25 degC (V):', pv.Vt())  # default is 25degC
print('Thermal voltage 300 K (V):', pv.Vt(300))
print('Silicon ni at 25 degC {:.3e} cm-3'.format(semi.ni_misiakos()))
print('Silicon ni at 300 K {:.3e} cm-3'.format(semi.ni_misiakos(300)))

print('n-type cell with doping level of 2e15')
ni = semi.ni_misiakos()
n0, p0 = semi.equilibrium_carrier(2e15, ni)  # n-type dioping at 1e15
print('Majority n {:.3e}, Minority p: {:.3e}'.format(n0, p0))
dn = 1e15
dEc, dEv = semi.bandgap_schenk(n0 + dn, p0 + dn, n0, p0, dn)
print('BGN Ec: {:.2e} eV, Ev {:.2e} eV'.format(dEc, dEv))

print('nieff {:.3e}'.format(semi.n_ieff(n0, p0, 1e15)))


