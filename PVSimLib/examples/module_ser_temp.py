from photovoltaic import *
import matplotlib.pyplot as plt
import numpy as np

# Series resistance losses
Rseries = 0.61  # ohm cm2
J_SC = 36.43
V_OC = 0.67
J_0 = 1.78e-10
FF = 0.84

print('************************* Station Info *************************')
city, station, GMT_offset, latitude, longitude, altitude, ETR, GHI, DNI, DHI, ambient_temperature = \
    common.station_info('722780TYA.csv')
print('station number: ', station)
print('City: ', city)
print('GMT: ', GMT_offset)
print('Latitude (degrees):', latitude)
print('Longitude (degrees):', longitude)
print('Altitude (m): ', altitude)

print('\n')
P_yearly = sum(DNI + DHI)

efficiency = solarcells.efficiency(V_OC, J_SC/1000, FF, 1)

J_L = systems.light_intensity(J_SC, DNI, DHI)

P_module, P_yearly_module = systems.power_module(J_L / 1000, J_0 / 1000, 0)
P_module_np = np.array(P_module)
efficiency_no_effects = systems.efficiency(P_module, DNI, DHI)

P_module_ser, P_yearly_module_ser = systems.power_module(J_L / 1000, J_0 / 1000, Rseries)
P_module_ser_np = np.array(P_module_ser)
efficiency_ser = systems.efficiency(P_module_ser_np, DNI, DHI)

# Temperature losses
#  temperature of the module
R_thermal = 0.041 # (K/W)
module_temperature_change, module_temperature = systems.module_temp(DNI, DHI, R_thermal, efficiency, ambient_temperature)
P_module_temp, P_yearly_module_temp = systems.power_module(J_L / 1000, J_0 / 1000, Rseries, module_temperature)

P_module_temp_np = np.array(P_module_temp)
efficiency_temp = systems.efficiency(P_module_temp_np, DNI, DHI)

print('Yearly Energy : {:.2f} kWh/m²'.format(P_yearly_module))
print('Yearly Energy  with Series Resistance: {:.2f} kWh/m²'.format(P_yearly_module_ser))
print('Yearly Energy (kWh/m²) with Series Resistance and Temp effect: {:.2f} kWh/m²'.format(P_yearly_module_temp))
print('Efficiency (%): {:.2f}'.format(1e5 * P_yearly_module/P_yearly))
print('Efficiency after Rseries (%): {:.2f}'.format(1e5 * P_yearly_module_ser/P_yearly))
print('Efficiency after Rseries and Temp (%): {:.2f}'.format(1e5 * P_yearly_module_temp/P_yearly))

# This prints the energy from the module
if 1:  # 1 - plot the data, 0 - turn plot off
    label = 'Power from Module'
    plt.figure(label)
    plt.title('Power from the Module')
    plt.plot(P_module[:200], label="No effects")
    plt.plot(P_module_ser[:200], label="Series")
    plt.plot(P_module_temp[:200], label="Temp + Series")
    plt.legend(loc='center right')
    plt.xlabel('hour of the year')
    plt.ylabel('Power in hour interval (kWh/m²)')
    plt.title('Power from the Modules')
    plt.savefig('module_power.png')
    plt.show()

if 1:  # 1 - plot the data, 0 - turn plot off
    plt.title('Yearly Efficiency')
    plt.plot(efficiency_no_effects[:200], label="No effects")
    plt.plot(efficiency_ser[:200], label="R Series")
    plt.plot(efficiency_temp[:200], label="Temp +  R Series")
    plt.legend(loc='center right')
    plt.xlabel('hour of the year')
    plt.ylabel('Efficiency(%)')
    plt.title('Yearly Efficiency')
    plt.savefig('yearly_efficiency.png')
    plt.show()
