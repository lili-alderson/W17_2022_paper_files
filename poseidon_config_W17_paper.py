# POSEIDON configuraton settings

import numpy as np

#***** Define physical constants *****#

R_J = 7.1492e7     # Radius of Jupiter (m)
M_J = 1.898e27     # Mass of Jupiter (kg)
R_E = 6.371e6      # Radius of Earth (m)
R_Sun = 6.957e8   # Radius of Sun (m)

#***** System parameters *****#

planet_name = 'WASP17b'

# Stellar properties
Band = 'K'            # Spectral band stellar magnitude measured in (options: 'J', 'H', 'K')
App_mag = 10.2       # Apparant magnitude of the star
R_s = 1.49*R_Sun     # Radius of star (m)
T_s = 6550.0          # Stellar effective temperature (K)
err_T_s = 100.0       # Error on known stellar effective temperature (K)
Met_s = -0.25           # Stellar metallicity [log10(Fe/H_star / Fe/H_solar)]
log_g_s = 4.2        # Stellar log surface gravity (cgs by convention)

# Planet properties
R_p = 1.87*R_J     # Radius of planet (m)
M_p = 0.78*M_J    # Mass of planet (kg)
g_0 = 3.948         # Gravitational field of planet: g = GM/R_p^2 (m/s^2)
a_p = 7.025*R_s      # Semi major axis of planetary orbit (m)
b_p = 0.36*R_s      # Impact parameter of planet orbit (m)
T_eq = 1700.0       # Equilibrium temperature (currently unused)

Trans_dur = 4.392    # Transit duration (hours)

#***** Atmosphere Setup *****#

N_D = 100           # Number of depths (layer centres) in atmosphere
P_max = 1.0e2       # Pressure at lowest altitude considered (bar)
P_min = 1.0e-7      # Pressure at highest altitude considered (bar)

#***** Wavelength grid *****#
 
wl_grid = 'constant R'  # Options: 'uniform' / 'constant R' / 'line-by-line'
wl_min = 0.28            # Minimum wavelength (um) (including errorbars!)
wl_max = 5.2            # Maximum wavelength (um) (including errorbars!)
N_wl = 1000             # Number of wavelength points for evalauting spectrum (uniform grid) - value here ignored if uniform not chosen
R = 4000                # Spectral resolution for evalauting spectrum (constant R)      

#***** Enable or disable program features *****#

mode = 'plot'   # forward_model / retrieve / plot

if (mode == 'forward_model'):
    
    single_model = True
    load_observations = True
    do_retrieval = False
    sim_retrieval = False
    produce_sim_data = False
    check_best_fit = False
    plot_retrieved_spectrum = False
    analyse_Bayes = False 
    make_corner = False
    
    skip_preload = False   # Ignore initialisation step (only for debugging!)

elif (mode == 'retrieve'):
    
    single_model = False
    load_observations = True
    do_retrieval = True
    sim_retrieval = False
    produce_sim_data = False
    check_best_fit = False
    plot_retrieved_spectrum = False
    analyse_Bayes = False
    make_corner = True
    
    skip_preload = False   # Ignore initialisation step (only for debugging!)
    
elif (mode == 'plot'):
    
    single_model = False
    load_observations = True
    do_retrieval = False
    sim_retrieval = False
    produce_sim_data = False
    check_best_fit = False
    plot_retrieved_spectrum = True
    analyse_Bayes = False
    make_corner = False
    
    skip_preload = False   # Ignore initialisation step (only for debugging!)

#***** Model settings *****#

PT_profile = 'isotherm'           # Options: isotherm / gradient / Madhu
cloud_model = 'MacMad17'        # Options: cloud-free / MacMad17 / Iceberg
cloud_type = 'deck_haze'               # Options: deck_haze / haze / deck (only applies if cloud_model != cloud-free)
cloud_dim = '1'                   # Options: 1 (uniform) / 2 (Patchy terminator)
stellar_contam = 'No'             # Options: No / one-spot

# Now specify which chemical species to include in this model
param_species = ['Na', 'K', 'H2O', 'CH4', 'CO', 'CO2']    # Chemical species with parametrised mixing ratios

# THESE SETTINGS FOR FORWARD MODELS
#***** Single model run parameters *****#

# (1) Mixing ratios (corresponding to included chemical species)

log_X_set = np.array([[-6.0, -7.0, -4.0]])     # Log abundances

# (2) P-T profile
    
# Specify P-T profile parameters
if (PT_profile == 'isotherm'):  
    PT_set = [1700]                  # T
elif (PT_profile == 'gradient'):  
    PT_set = [800, 2000]            # T_high, T_deep
elif (PT_profile == 'Madhu'):     
    PT_set = [1.2, 1.0, -2.0, -5.0, 0.0, 2200.0]  # a1, a2, log(P1,2,3), T_deep

# Radius at reference pressure
R_p_ref_set = 0.96*(R_p/R_J)  

# (3) Cloud properties
    
# 1D model (uniform clouds)
if (cloud_dim == '1'):  
    if (cloud_model == 'cloud-free'):   # No cloud parameters for clear atmosphere
        clouds_set = []
    else:            
        if (cloud_model == 'MacMad17'):
            if (cloud_type == 'deck'):      
                clouds_set = [-2.0]               # log(P_cloud)
            elif (cloud_type == 'haze'):      
                clouds_set = [2.0, -8.0]          # log(a), gamma
            elif (cloud_type == 'deck_haze'):
                clouds_set = [2.0, -7.0, -2.0]    # log(a), gamma, log(P_cloud)
            
# 2D model (patchy clouds)
elif (cloud_dim == '2'):  
    if (cloud_model == 'cloud-free'):   # No cloud parameters for clear atmosphere
        clouds_set = []
    else:
        if (cloud_model == 'MacMad17'):
            if (cloud_type == 'deck'):      
                clouds_set = [-2.0, 0.5, 45.0]               # log(P_cloud), phi_c, phi_0
            elif (cloud_type == 'deck_haze'):
                clouds_set = [2.0, -7.0, -2.0, 0.5, 45.0]    # log(a), gamma, log(P_cloud), phi_c, phi_0
            elif (cloud_type == 'haze'):      
                clouds_set = [2.0, -8.0, 0.5, 45.0]          # log(a), gamma, phi_c, phi_0

# (4) Stellar contamination parameters
            
if (stellar_contam == 'No'):    # No stellar contamination
    stellar_set = []
elif (stellar_contam == 'one-spot'):
    stellar_set = [0.05, 6800.0, T_s]     # f_het, T_het, T_phot
    
    
#***** Data sources and settings *****#

# Specify file locations of each data file
STIS_G430 = planet_name + '_STIS_G430.dat'
STIS_G750 = planet_name + '_STIS_G750.dat'
WFC3_G102 = planet_name + '_WFC3_G102.dat'
WFC3_G141 = planet_name + '_WFC3_G141.dat'
IRAC1 = planet_name + '_Spitzer_IRAC1.dat'
IRAC2 = planet_name + '_Spitzer_IRAC2.dat'

# Provide full list of instruments to use
instruments = np.array(['STIS_430', 'STIS_G750', 'WFC3_G102', 'WFC3_G141', 'IRAC1', 'IRAC2'])   

# Specify datasets above to use
datasets = [STIS_G430, STIS_G750, WFC3_G102, WFC3_G141, IRAC1, IRAC2]    


#***** Retrieval settings *****#

# Set lower prior limits for parameters
prior_lower_X = np.array([-12.0, -3.0])                # Lower prior for log(X_i) and delta log(X_i)
prior_lower_R_p_ref = 0.85*(R_p/R_J)                   # Lower prior for R_p_ref (in Jupiter radii)
prior_lower_stellar = np.array([0.0, 0.6*T_s])         # Lower priors for f, T_het

# Set upper prior limits for parameters
prior_upper_X = np.array([-1.0, 3.0])                  # Upper prior for log(X_i) and delta log(X_i)
prior_upper_R_p_ref = 1.15*(R_p/R_J)                   # Upper prior for R_p_ref (in Jupiter radii)
prior_upper_stellar = np.array([0.5, 1.2*T_s])         # Upper priors for f, T_het

# Set mean and std for Gaussian parameter priors
prior_gauss_T_phot = np.array([T_s, err_T_s])          # Gaussian prior for T_phot

# Cloud parameter priors
if (cloud_model == 'MacMad17'):
    prior_lower_clouds = np.array([-4.0, -20.0, -6.0, 0.0, -180.0])     # Lower priors for log(a), gamma, log(P_cloud), phi_c, phi_0
    prior_upper_clouds = np.array([8.0, 2.0, 2.0, 1.0, 180.0])          # Upper priors for log(a), gamma, log(P_cloud), phi_c, phi_0
else:
    prior_lower_clouds = np.array([])
    prior_upper_clouds = np.array([])

# P-T profile priors
T_min = 400.0      # Minimum temperature to be considered
T_max = 2300.0     # Maximum temperature to be considered

if (PT_profile == 'isotherm'):
    prior_lower_PT = np.array([T_min])      # Lower prior for T
    prior_upper_PT = np.array([T_max])      # Upper prior for T
elif (PT_profile == 'gradient'):
    prior_lower_PT = np.array([T_min, T_min, 0.0])        # Lower priors for T_high, T_deep, delta_T
    prior_upper_PT = np.array([T_max, T_max, 1000.0])     # Upper priors for T_high, T_deep, delta_T    
elif (PT_profile == 'Madhu'):
    prior_lower_PT = np.array([0.02, 0.02, -6.0, -6.0, -2.0, T_min])   # Lower priors for a1, a2, log(P1,2,3), T_deep
    prior_upper_PT = np.array([2.00, 2.00, 1.0, 1.0, 1.0, T_max])      # Upper priors for a1, a2, log(P1,2,3), T_deep


# Specify sampling algorithm and properties
sampling_algorithm = 'MultiNest'                          # Options: MultiNest
N_live = 4000                                             # Number of MultiNest live points
base_name = planet_name + '_' + str(N_live) + 'atmo_2EB-'    # Output file base name KEEP - ON THE END!!
# KEEP THE BASENAME SHORT!!!!!!!




#**************************************************************#
#***** Advanced settings (shouldn't need to change these) *****#
#**************************************************************#


spectrum_type = 'transmission'    # Options: transmission 
rad_transfer = 'geometric'        # Options: geometric

X_dim = '1'                       # Options: 1 (uniform) / 2 (Evening-Morning or Day-Night) / 3 (Evening-Morning-Day-Night)
PT_dim = '1'                      # Options: 1 (uniform) / 2 (Evening-Morning or Day-Night) / 3 (Evening-Morning-Day-Night)
X_profile = 'isochem'             # Options: isochem / gradient
TwoD_type = 'E-M'                 # Options: E-M (2D Evening-Morning model) / D-N (2D Day-Night model)
TwoD_param_scheme = 'difference'  # Options: absolute / difference
term_transition = 'linear'        # Options: linear (linear terminator transition - 2D / 3D models only)
He_fraction_setting = 'fixed'     # Options: fixed / free (if set to free, then X_He becomes a free parameter)
chemistry_prior = 'log-uniform'   # Options: log-uniform (Bulk component known) / CLR (a priori unknown, max 8 gases)
offsets_applied = 'No'            # Options: No / relative
error_inflation = 'No'            # Options: No / Line_2015

# Specify the chemical species for which we consider non-uniform gradients 
species_EM_gradient = []
species_DN_gradient = []
species_vert_gradient = []

# Species filling most of atmosphere
bulk_species = ['H2', 'He']           

# If fixing He/H2 ratio, use this value
He_fraction = 0.17

# Other pressure settings
P_deep = 10     # 'Anchor' pressure below which the atmosphere is homogenous
P_high = 1.0e-5   # Pressure where temperature parameters are defined ('gradient' PT profile)

# Geometry settings
alpha_set = 0.1   # Angular width of Evening-Morning terminator transition (2D E-M and 3D only)
beta_set = 0.1    # Angular width of Day-Night terminator transition (2D D-N and 3D only)

# Resolution for integration along / across terminator (if term_transition = 'linear')
N_slice_EM = 2    # Number of azimuthal slices across Evening-Morning terminator transition (Even number !)
N_slice_DN = 2    # Number of angular slices along Day-Night terminator transition(Even number !)

# Other settings
offsets_set = []      
offset_datasets = []  
err_inflation_set = []    

# Absorption.py settings
opacity_treatment = 'Opacity-sample'   # Options: 'line-by-line' / 'Opacity-sample'
opacity_database = 'High-T'            # Options: 'High-T' / 'Temperate'

line_by_line_resolution = 0.01    # If using line-by-line mode, specify wavenumber resolution of cross sections calculation

T_fine_step = 10.0      # Temperature resolution for pre-interpolation of opacities (K)
T_fine_min = T_min      # Minimum temperature on fine temperature grid (K)
T_fine_max = T_max      # Maximum temperature on fine temperature grid (K)

N_D_pre_inp = 40        # Pressure resolution for pre-interpolation of opacities

# Stellar.py settings
T_phot_min = T_s - 10.0*err_T_s     # Minimum T_phot on pre-computed grid (-10 sigma)
T_phot_max = T_s + 10.0*err_T_s     # Maximum T_phot on pre-computed grid (+10 sigma)
T_phot_step = err_T_s/10.0          # T_phot pre-computed grid resolution (0.1 sigma)
T_het_min = 0.6*T_s                 # Minimum T_het on pre-computed grid
T_het_max = 1.2*T_s                 # Maximum T_het on pre-computed grid
T_het_step = 10.0                   # T_het pre-computed grid resolution (K)

# Other MultiNest settings
sampling_target = 'parameter'       # Options: parameter / model 
ev_tol = 0.5                        # Evidence tolerance factor

prior_lower_geometry = np.array([0.1, 0.1])            # Lower prior for alpha and beta (degrees)
prior_lower_offsets = np.array([-1.0e-3, -1.0e-3])     # Lower priors for linear DC offsets and relative offset 

prior_upper_geometry = np.array([180.0, 60.0])         # Upper prior for alpha and beta (degrees)
prior_upper_offsets = np.array([1.0e-3, 1.0e-3])       # Upper priors for linear DC offsets and relative offset 

#***** Synthetic data generation *****#

# Settings for simple constant precision constant resolution run
std_data = 100       # Standard deviation of synthetic data (ppm)
R_data = 40          # Spectral resolution of synthetic data (wl/delta_wl)
wl_data_min = 0.44   # Starting wavelength of synthetic data range (um)
wl_data_max = 0.54   # Ending wavelength of synthetic data range (um)

#std_data = 100       # Standard deviation of synthetic data (ppm)
#R_data = 40          # Spectral resolution of synthetic data (wl/delta_wl)
#wl_data_min = 0.58   # Starting wavelength of synthetic data range (um)
#wl_data_max = 0.97   # Ending wavelength of synthetic data range (um)

#std_data = 50        # Standard deviation of synthetic data (ppm)
#R_data = 60          # Spectral resolution of synthetic data (wl/delta_wl)
#wl_data_min = 1.10   # Starting wavelength of synthetic data range (um)
#wl_data_max = 1.70   # Ending wavelength of synthetic data range (um)