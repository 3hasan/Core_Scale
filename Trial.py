                                            ### Core Scale Model for Deep Bed Filtration ###
                                            

import numpy as np
import matplotlib.pyplot as plt

### SOME CONSTANTS

g = 980 # cm/min2
particle_density = 2.65 # g/cm3
porosity = 0.25 
permeability = 2 #Darcy
reaction_rate = 3

### FLUID PROPERTIES

fluid_density = 1 # g/cm3
fluid_dynamic_viscosity = 0.89 #cp

###CORE DIMENSIONS

Core_Length = 5        # cm
Core_Area = 5.23       # cm2
Core_Height = 2 * np.sqrt(Core_Area / np.pi) # cm
Core_Volume = Core_Length * Core_Area #cm3

### INJECTION PROPERTIES

injection_total_concentration = 20  # ppm
injection_rate = 5 # cm3/min
darcy_velocity = injection_rate / Core_Area  # cm/min

### TIME DISCRETIZATION

time_total = 60 #minutes   
time_step = 0.01  #minutes          
N_time_steps = int(time_total / time_step)
time = np.linspace(0, time_total, N_time_steps)

### SPACE DISCRETIZATION

N_grid_cells = 100
deltaX = Core_Length / N_grid_cells # cm

Cell_Volume = deltaX * Core_Area #cm3
cell_centers = np.linspace(-deltaX / 2, Core_Length + deltaX / 2, N_grid_cells + 2)

concentration = np.zeros((N_grid_cells + 2, N_time_steps))
concentration[0, :] = injection_total_concentration # Boundary cond.

def update_concentration(concentration, darcy_velocity, deltaX, reaction_rate, time_step, N_grid_cells, j):
    concentration[-1, j] = concentration[-2, j-1] # Boundary cond.

    for i in range(1, N_grid_cells + 1): #Central difference
        concentration[i, j] = concentration[i, j-1] + time_step * (-darcy_velocity * (concentration[i+1, j-1] - concentration[i-1, j-1]) / (2 * deltaX) - reaction_rate * concentration[i, j-1])

        
for j in range(1, N_time_steps):
    update_concentration(concentration, darcy_velocity, deltaX, reaction_rate, time_step, N_grid_cells, j)


plt.scatter(cell_centers, concentration[:, -1])
plt.ylabel('Concentration')
plt.xlabel('X centimeters')
plt.show()



 

























































































                ### Data from 'BCR Information - Certification report on reference materials of defined particle size', (1980), Commission of the European Communities ###

#particle_diameters = np.array([20, 11.3, 9.5, 8, 6.7, 5.7, 4.8, 3.5, 2.9, 2.4, 2, 1.7, 1.4, 1.2]) * 1e-3  # microns to cm
#cumulative_percentages = np.array([0.998, 0.976, 0.954, 0.924, 0.881, 0.832, 0.768, 0.602, 0.495, 0.374, 0.274, 0.22, 0.134, 0.10]) # 0 to 1

#                           ### DATA MANUALLY EXTRACTED FROM 'Bentheimer sandstone revisited for experimental purposes' Peksa et al.(2015) ###


### Pore Throat Sizes(milimeters) and Cumulative Distributions (%)

#throat_sizes = np.array([0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.028, 0.032, 0.036, 0.040, 0.044, 0.048, 0.052, 
#                         0.056, 0.060, 0.064, 0.068, 0.072, 0.076, 0.080, 0.084, 0.088, 0.092, 0.096, 0.100, 0.104, 0.108])

#throat_size_cumulative_percentage = np.array([99.09, 96.65, 85.87, 65.22, 49.64, 38.41, 30.53, 25.00, 
#                                              20.20, 16.58, 13.68, 10.96, 9.15, 7.34, 5.71, 4.98, 3.89, 2.99, 2.45, 2.08, 1.54, 1.36, 1.00, 0.72, 0.45, 0.27, 0.25])

### Pore Body Sizes(milimeters) and Cumulative Distributions (%)

#body_sizes = np.array([0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.028, 0.032, 0.036, 0.040, 0.044, 
#                       0.048, 0.052, 0.056, 0.060, 0.064, 0.068, 0.072, 0.076, 0.080, 0.084, 0.088, 0.092, 0.096, 0.100, 0.104, 0.108, 0.112, 0.116, 0.120, 
#                       0.124, 0.128, 0.132, 0.136, 0.140, 0.144, 0.148, 0.152, 0.156, 0.160, 0.164, 0.168, 0.172, 0.176, 0.180, 0.184, 0.188, 0.192, 0.196, 0.200, 0.204, 
#                       0.208, 0.212, 0.216, 0.220, 0.224, 0.228, 0.232, 0.236, 0.240, 0.244, 0.248, 0.252, 0.256, 0.260])

#body_size_cumulative_percentage = np.array([99.90, 99.82, 87.86, 59.60, 41.12, 30.43, 24.73, 21.56, 
#                                            20.11, 19.11, 18.30, 17.75, 17.21, 16.85, 16.67, 16.12, 15.76, 15.49, 15.04, 14.76, 14.40, 14.22, 14.13, 
#                                            13.41, 13.13, 12.77, 12.50, 12.23, 11.96, 11.78, 11.23, 10.96, 10.60, 10.05, 9.87, 9.33, 8.97, 8.70, 8.24, 7.79, 7.70, 6.97, 
#                                            6.70, 6.16, 5.89, 5.62, 5.34, 5.07, 4.80, 4.44, 3.99, 3.71, 3.62, 3.35, 3.17, 3.08, 2.81, 2.54, 2.08, 1.99, 1.81, 1.54, 1.45, 
#                                            1.18, 0.63])


### Grain Sizes(milimeters) and Cumulative Distributions (%)

#grain_sizes = np.array([0.0443, 0.0567, 0.0691, 0.0814, 0.0938, 0.1061, 0.1185, 0.1308, 0.1432, 0.1556, 
#                        0.1679, 0.1803, 0.1926, 0.2050, 0.2173, 0.2297, 0.2421, 0.2544, 0.2668, 0.2791, 0.2915, 0.3038, 0.3162, 0.3286, 0.3409, 
#                        0.3533, 0.3656, 0.3780, 0.3903, 0.4027, 0.4151, 0.4274, 0.4398, 0.4521, 0.4645, 0.4768, 0.4892, 0.5015, 0.5139, 0.5263, 
#                        0.5386, 0.5510, 0.5633, 0.5757, 0.5880, 0.6004, 0.6128, 0.6251])

#grain_size_cumulative_percentage = np.array([99.99, 99.95, 99.85, 99.78, 99.68, 99.53, 99.28, 98.89, 
#                                             98.31, 97.52, 96.46, 95.05, 93.06, 90.20, 86.21, 80.94, 74.34, 67.11, 59.96, 53.33, 47.27, 41.75, 36.77, 
#                                             32.30, 28.33, 24.84, 21.82, 19.24, 17.10, 15.34, 13.91, 12.77, 11.88, 11.18, 10.62, 10.16, 9.75, 9.35, 
#                                             8.90, 8.36, 7.67, 6.80, 5.75, 4.62, 3.54, 2.63, 2.01, 1.81])


##               ### REACTION PROCESSES ###

## Sedimentation
#def stokes_settling_velocity(particle_diameters, particle_density, fluid_density, fluid_viscosity):
#    v_settling = ((particle_diameters**2) * g * (particle_density - fluid_density)) / (18 * (fluid_dynamic_viscosity / fluid_density))
#    return v_settling 

## Interception

## Direct Impaction