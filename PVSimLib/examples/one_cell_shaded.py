from photovoltaic import *
import matplotlib.pyplot as plt

T = 300
Rseries = 0.61  # ohm cm2
J_SC = 36.43
V_OC = 0.67
J_0 = 1.78e-10
FF = 0.84

efficiency = solarcells.efficiency(V_OC, J_SC/1000, FF, 1)
shading = 0.4
number_series = 8  # Number UNSHADED in series
shading_conc = (number_series + shading)/(number_series + 1)
number_points = 100  # Number of points in each IV curve
Vbr = -0.5

city, station, GMT_offset, latitude, longitude, altitude, ETR, GHI, DNI, DHI, ambient_temperature = \
    common.station_info('722780TYA.csv')
print('station number: ', station)
print('City: ', city)
print('GMT: ', GMT_offset)
print('Latitude (degrees):', latitude)
print('Longitude (degrees):', longitude)
print('Altitude (m): ', altitude)

R_thermal = 0.041 # (K/W)
module_temperature_change, module_temperature = systems.module_temp(DNI, DHI, R_thermal, efficiency, ambient_temperature)

Voc_total, Jsc_total, FF_total, Vmp_total, Jmp_total, Pmp_total, voltage_total, current = \
    systems.IV_curve_module(J_SC / 1000, J_0 / 1000, Rseries, shading, Vbr, number_series, T)

if 1:  # 1 - plot the data, 0 - turn plot off
    label = 'Module Shaded Series IV Curves '
    plt.title('IV Curves for Shaded Series Cells with {:.0f} in series'.format(number_series + 1))
    plt.plot(voltage_total, current, label='Shading {:.0f} %'.format(shading*100))
    plt.scatter(Vmp_total, Jmp_total)
    plt.text(Vmp_total,Jmp_total, ' Maximum \nPower Point')
    plt.xlim(-0.5)
    plt.xlabel('Voltage (V)')
    plt.ylabel('Current Density (A/cm²)')
    plt.legend(loc='center left')
    plt.ylim(0)
    plt.twinx()
    plt.plot(voltage_total, voltage_total*current, c='red')
    plt.ylim(0)
    plt.ylabel('Power Density(W/cm²)')
    plt.show()

J_L = systems.light_intensity(J_SC, DNI, DHI)
P_yearly_module_unshaded, P_yearly_module_shaded = systems.power_module_shaded(J_L / 1000, J_0 / 1000, Rseries,
                                                                               shading, Vbr, number_series, T)
P_yearly = sum(DNI+DHI)
print('\n')
print('Energy from Sun over year', P_yearly / 1000)
print('Energy from Module {:.2f} kWh/m²'.format(P_yearly_module_unshaded))
print('Energy from Module with shading {:.2f} kWh/m²'.format(P_yearly_module_shaded))
print('Average Efficiency = {:.2f}% '.format(1e5 * P_yearly_module_shaded / P_yearly))
print('\n')

shading = 0.1
P_yearly_shaded = []
shading_vec = []

while(shading < 0.6):
    P_yearly_module_unshaded, P_yearly_module_shaded = systems.power_module_shaded(J_L / 1000, J_0 / 1000, Rseries,
                                                                                   shading, Vbr, number_series, T)
    P_yearly_shaded.append(P_yearly_module_shaded)
    shading_vec.append(shading*100)
    shading = shading + 0.01

if 1:  # 1 - plot the data, 0 - turn plot off
    plt.title('IV Curves for Shaded Series Cells with {:.0f} in series'.format(number_series + 1))
    plt.plot(shading_vec,P_yearly_shaded)
    plt.grid(True)
    plt.xlim(5)
    plt.xlabel('Shading (%)')
    plt.ylabel('Energy from Module with shading (kWh/m)')
    plt.show()