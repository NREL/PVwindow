from photovoltaic import *

W = 15.6                    # Solar cell width (cm)
Z = 15.6                    # Solar cell height (cm)

# Emitter parameters
We = 0.25e-4        # emitter thickness (cm)

# Base Parameters
Wb = 250e-4         # Base thickness (cm)

# Concentration parameters
Nd = 2e19           # Majority carrier concentration on N side
Na = 3e16           # Majority carrier concentration on P side

# Front contacts parameters
bus_resistivity = 1.6e-6    # Bus resistivity (Ohm*cm)
bus_width = 0.2             # Bus width (cm)
bus_height = 400e-4         # Bus height (cm)
bus_number = 3              # Number of buses

finger_resistivity = 4.8e-6                 # Finger resistivity (Ohm*cm)
finger_spacing = 0.2                        # Finger spacing (cm)
finger_width = 120e-4                       # Finger width (cm)
finger_length = Z/(2*bus_number)            # Finger length (cm)
finger_height = 20e-4                       # Finger height (cm)
finger_number = round(W/finger_spacing)     # Number of fingers

# Series Resistance parameters
resistivity_base = semiconductors.resistivity_Si_p(Na)                         # Base resistivity (Ohm*cm)
resistivity_emitter = semiconductors.resistivity_Si_n(Nd)                      # Emitter resistivity (Ohm*cm)
R_sheet = solarcells.sheet_resistance(resistivity_emitter, We)    # Emitter sheet resistivity (Ohm/square)


Jmp = 33.41
Vmp = 0.62

print('************************* Losses *************************')
print('Base resistivity: {:.2f} Ohm*cm'.format(resistivity_base))
base_resistance = solarcells.base_resistance(resistivity_base, Wb)
print('Base Resistance: {:.2f} Ohm*cm^2'.format(base_resistance))
base_losses = solarcells.base_losses(resistivity_base, Wb, Vmp, Jmp)
print('Base Losses: {:.2f} %'.format(base_losses))

print('\n')
print('Emitter resistivity: {:.2e} Ohm*cm'.format(resistivity_emitter))
print('Sheet resistance: {:.2f} Ohm/square'.format(R_sheet))
emitter_resistance = solarcells.emitter_resistance(R_sheet, finger_spacing)
print('Emitter Resistance: {:.2f} Ohm*cm^2'.format(emitter_resistance))
emitter_losses = solarcells.finger_sheet(finger_spacing, Jmp/1000, R_sheet, Vmp)
print('Emitter Losses: {:.2f} %'.format(emitter_losses))

print('\n')
finger_resistance = solarcells.finger_resistance(finger_resistivity, finger_spacing, finger_width, finger_length,
                                                 finger_height)
print('Finger Resistance: {:.2f} Ohm*cm^2'.format(finger_resistance))
finger_losses = solarcells.finger_resistivity(finger_length, Jmp*0.001, finger_spacing, finger_resistivity,
                                              finger_width, finger_height, Vmp)
print('Finger Losses: {:.2f} %'.format(finger_losses))

print('\n')
bus_resistance = solarcells.bus_resistance(W, bus_resistivity, bus_number, bus_width, bus_height)
print('Bus Resistance: {:.2f} Ohm*cm^2'.format(bus_resistance))
bus_losses = solarcells.bus_losses(W, Z, bus_resistivity, Jmp*0.001, Vmp, bus_number, bus_width, bus_height)
print('Bus Losses: {:.2f} %'.format(bus_losses))

print('\n')
transparency = solarcells.transparency(W, Z, bus_number, finger_number, finger_width, bus_width)
print('transparency: {:.2f} %'.format(transparency))


print('\n')
print('************************* Series Resistance *************************')
Rs = bus_resistance + finger_resistance + base_resistance + emitter_resistance
print('Rseries: {:.2f} Ohm cm^2'.format(Rs))
