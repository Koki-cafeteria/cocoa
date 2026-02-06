import sys, platform, os
os.environ['OMP_NUM_THREADS'] = '12'
import matplotlib
import math
from matplotlib import pyplot as plt
import numpy as np
#import euclidemu2
import scipy

sys.path.append("/home/tanida/cocoa_v41/cocoa/Cocoa/projects/mg_hsc_y3/interface/")
import cosmolike_mg_hsc_y3_interface as ci
from getdist import IniFile
import itertools
import iminuit
import functools
print(sys.version)
print(os.getcwd())

# GENERAL PLOT OPTIONS
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['xtick.bottom'] = True
matplotlib.rcParams['xtick.top'] = False
matplotlib.rcParams['ytick.right'] = False
matplotlib.rcParams['axes.edgecolor'] = 'black'
matplotlib.rcParams['axes.linewidth'] = '1.0'
matplotlib.rcParams['axes.labelsize'] = 'medium'
matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['grid.linewidth'] = '0.0'
matplotlib.rcParams['grid.alpha'] = '0.18'
matplotlib.rcParams['grid.color'] = 'lightgray'
matplotlib.rcParams['legend.labelspacing'] = 0.77
matplotlib.rcParams['savefig.bbox'] = 'tight'
matplotlib.rcParams['savefig.format'] = 'pdf'
matplotlib.rcParams['text.usetex'] = True

# Jupyter Notebook Display options
import IPython
IPython.display.display(IPython.display.HTML("<style>:root { --jp-notebook-max-width: 85% !important; }</style>"))
IPython.display.display(IPython.display.HTML("<style>div.output_scroll { height: 54em; }</style>"))




from scipy.interpolate import InterpolatedUnivariateSpline as iuspline

# IMPORT CAMB
sys.path.insert(0, "/home/tanida/cocoa_v41/cocoa/Cocoa"+"/external_modules/code/MGCAMB/build/lib")
import camb
from camb import model
print('Using CAMB %s installed at %s'%(camb.__version__,os.path.dirname(camb.__file__)))


CAMBAccuracyBoost = 1.1
non_linear_emul = 2
CLprobe="xi"

path= "/home/tanida/cocoa_v41/cocoa/Cocoa/projects/mg_hsc_y3/data/"
data_file="mg_hsc_y3_M1_GGL0.05.dataset"
# path= "/home/tanida/cocoa_v41/cocoa/Cocoa/projects/lsst_y1/data/"
# data_file="lsst_y1_M1_GGL0.05.dataset"

IA_model = 0
IA_redshift_evolution = 0


mu0 = 0.0
sigma0 =0.0
As_1e9 = 2.2065
ns = 0.9645
omegam = 1.- 0.6844
H0 = np.sqrt((0.02225+0.1198)/omegam)*100.
omegab = 0.02225/((0.02225+0.1198)/omegam)
mnu = 0.06
HSC_DZ_S1 = 0.0
HSC_DZ_S2 = 0.0
HSC_DZ_S3 = 0.0
HSC_DZ_S4 = 0.0
HSC_DZ_S5 = 0.0
HSC_M1 = 0.0
HSC_M2 = 0.0
HSC_M3 = -0.0
HSC_M4 = -0.0
HSC_M5 = -0.0
HSC_A1_1 = 0.0
HSC_A1_2 = 0.0
w0pwa = -1.0
w = -1.0


def get_camb_cosmology(omegam = omegam, omegab = omegab, H0 = H0, ns = ns, 
                       As_1e9 = As_1e9, w = w, w0pwa = w0pwa, AccuracyBoost = 1.0, 
                       kmax = 10, k_per_logint = 20, CAMBAccuracyBoost=1.1):

    As = lambda As_1e9: 1e-9 * As_1e9
    wa = lambda w0pwa, w: w0pwa - w
    omegabh2 = lambda omegab, H0: omegab*(H0/100)**2
    omegach2 = lambda omegam, omegab, mnu, H0: (omegam-omegab)*(H0/100)**2-(mnu*(3.046/3)**0.75)/94.0708
    omegamh2 = lambda omegam, H0: omegam*(H0/100)**2

    CAMBAccuracyBoost = CAMBAccuracyBoost*AccuracyBoost
    kmax = max(kmax/2.0, kmax*(1.0 + 3*(AccuracyBoost-1)))
    k_per_logint = max(k_per_logint/2.0, int(k_per_logint) + int(3*(AccuracyBoost-1)))
    extrap_kmax = max(max(2.5e2, 3*kmax), max(2.5e2, 3*kmax) * AccuracyBoost)

    z_interp_1D = np.concatenate( (np.concatenate( (np.linspace(0,2.0,1000),
                                                    np.linspace(2.0,10.1,200)),
                                                    axis=0
                                                 ),
                                   np.linspace(1080,2000,20)),
                                   axis=0)
    
    z_interp_2D = np.concatenate((np.linspace(0, 2.0, 95), np.linspace(2.25, 10, 5)),  axis=0)

    log10k_interp_2D = np.linspace(-4.2, 2.0, 1200)

    pars = camb.set_params(H0=H0, 
                           ombh2=omegabh2(omegab, H0), 
                           omch2=omegach2(omegam, omegab, mnu, H0), 
                           mnu=mnu, 
                           omk=0, 
                           tau=0.06,  
                           As=As(As_1e9), 
                           ns=ns, 
                           halofit_version='takahashi', 
                           lmax=10,
                           AccuracyBoost=CAMBAccuracyBoost,
                           lens_potential_accuracy=1.0,
                           num_massive_neutrinos=1,
                           nnu=3.046,
                           accurate_massive_neutrino_transfers=False,
                           k_per_logint=k_per_logint,
                           kmax = kmax);
    
    pars.set_dark_energy(w=w, wa=wa(w0pwa, w), dark_energy_model='ppf');    
    
    pars.NonLinear = model.NonLinear_both
    
    pars.set_matter_power(redshifts = z_interp_2D, kmax = kmax, silent = True);
    results = camb.get_results(pars)
    
    PKL  = results.get_matter_power_interpolator(var1="delta_tot", var2="delta_tot", nonlinear = False, 
                                                 extrap_kmax = extrap_kmax, hubble_units = False, k_hunit = False);
    
    PKNL = results.get_matter_power_interpolator(var1="delta_tot", var2="delta_tot",  nonlinear = True, 
                                                 extrap_kmax = extrap_kmax, hubble_units = False, k_hunit = False);
    
    lnPL = np.empty(len(log10k_interp_2D)*len(z_interp_2D))
    for i in range(len(z_interp_2D)):
        lnPL[i::len(z_interp_2D)] = np.log(PKL.P(z_interp_2D[i], np.power(10.0,log10k_interp_2D)))
    lnPL  += np.log(((H0/100.)**3)) 
    
    lnPNL  = np.empty(len(log10k_interp_2D)*len(z_interp_2D))
    if non_linear_emul == 1:
        params = { 'Omm'  : omegam, 
                   'As'   : As(As_1e9), 
                   'Omb'  : omegab,
                   'ns'   : ns, 
                   'h'    : H0/100., 
                   'mnu'  : mnu,  
                   'w'    : w, 
                   'wa'   : wa(w0pwa, w)
                 }
        kbt, bt = euclidemu2.get_boost( params, 
                                        z_interp_2D, 
                                        np.power(10.0, np.linspace( -2.0589, 0.973, len(log10k_interp_2D)))
                                      )
        log10k_interp_2D = log10k_interp_2D - np.log10(H0/100.)
        
        for i in range(len(z_interp_2D)):    
            lnbt = scipy.interpolate.interp1d(np.log10(kbt), np.log(bt[i]), kind = 'linear', 
                                              fill_value = 'extrapolate', 
                                              assume_sorted = True)(log10k_interp_2D)
            lnbt[np.power(10,log10k_interp_2D) < 8.73e-3] = 0.0
            lnPNL[i::len(z_interp_2D)]  = lnPL[i::len(z_interp_2D)] + lnbt
    elif non_linear_emul == 2:
        for i in range(len(z_interp_2D)):
            lnPNL[i::len(z_interp_2D)] = np.log(PKNL.P(z_interp_2D[i], np.power(10.0, log10k_interp_2D)))            
        log10k_interp_2D = log10k_interp_2D - np.log10(H0/100.)
        lnPNL += np.log(((H0/100.)**3))

    G_growth = np.sqrt(PKL.P(z_interp_2D,0.0005)/PKL.P(0,0.0005))
    G_growth = G_growth*(1 + z_interp_2D)/G_growth[len(G_growth)-1]

    chi = results.comoving_radial_distance(z_interp_1D, tol=1e-4) * (H0/100.)

    return (log10k_interp_2D, z_interp_2D, lnPL, lnPNL, G_growth, z_interp_1D, chi)


def DeltaSigma(nR = 30,
           R_min = 5.65421091e-02, 
           R_max = 7.07437353e+01,
           ntheta = 30,
           theta_min_arcmin = 2.822145388443246250e-01, 
           theta_max_arcmin = 3.189063198818608384e+02,
           omegam = omegam, 
           omegab = omegab, 
           H0 = H0, 
           ns = ns, 
           As_1e9 = As_1e9, 
           w = w, 
           w0pwa = w0pwa,
           #A1  = [HSC_A1_1, HSC_A1_2, 0, 0, 0], 
           #A2  = [0, 0, 0, 0, 0],
           #BTA = [0, 0, 0, 0, 0],
           #shear_photoz_bias = [HSC_DZ_S1, HSC_DZ_S2, HSC_DZ_S3, HSC_DZ_S4, HSC_DZ_S5],
           #M = [HSC_M1, HSC_M2, HSC_M3, HSC_M4, HSC_M5],
           baryon_sims = None,
           AccuracyBoost = 1.0, 
           kmax = 1000, 
           k_per_logint = 20, 
           CAMBAccuracyBoost=1.1,
           CLAccuracyBoost = 1.0, 
           CLIntegrationAccuracy = 1):

    (log10k_interp_2D, z_interp_2D, lnPL, lnPNL, G_growth, z_interp_1D, chi) = get_camb_cosmology(omegam=omegam, 
        omegab=omegab, H0=H0, ns=ns, As_1e9=As_1e9, w=w, w0pwa=w0pwa, AccuracyBoost=AccuracyBoost, kmax=kmax,
        k_per_logint=k_per_logint, CAMBAccuracyBoost=CAMBAccuracyBoost)

    CLAccuracyBoost = CLAccuracyBoost * AccuracyBoost
    CLSamplingBoost = CLAccuracyBoost * AccuracyBoost
    CLIntegrationAccuracy = max(0, CLIntegrationAccuracy + 5*(AccuracyBoost-1.0))
    ci.init_accuracy_boost(1.0, CLAccuracyBoost, int(CLIntegrationAccuracy))
    
    #ci.init_binning_flat(int(nR), R_min, R_max)
    ci.init_binning(int(ntheta), theta_min_arcmin, theta_max_arcmin)
    
    ci.set_cosmology(omegam = omegam, 
                     H0 = H0, 
                     log10k_2D = log10k_interp_2D, 
                     z_2D = z_interp_2D, 
                     lnP_linear = lnPL,
                     lnP_nonlinear = lnPNL,
                     G = G_growth,
                     z_1D = z_interp_1D,
                     chi = chi)
    #ci.set_nuisance_shear_calib(M = M)
    #ci.set_nuisance_shear_photoz(bias = shear_photoz_bias)
    #ci.set_nuisance_ia(A1 = A1, A2 = A2, B_TA = BTA)

    if baryon_sims is None:
        ci.reset_bary_struct()
    else:
        ci.init_baryons_contamination(sim = baryon_sims)
        
    DS = ci.DeltaSigma_flatsky() 
    #DS = ci.xi_pm_tomo()
    #DS = ci.w_gg_tomo()
    return (ci.get_binning_real_space_flat(), DS)

# Init Cosmolike
ini = IniFile(os.path.normpath(os.path.join(path, data_file)))

lens_file = ini.relativeFileName('nz_lens_file')

source_file = ini.relativeFileName('nz_source_file')


lens_ntomo = ini.int("lens_ntomo")

source_ntomo = ini.int("source_ntomo")
#source_ntomo = 3

ci.initial_setup()

ci.init_accuracy_boost(1.0, 1.0, int(1))

ci.init_cosmo_runmode(is_linear = False)

ci.init_source_sample(filename = source_file, 
                       ntomo_bins = int(source_ntomo))
                      #ntomo_bins = 3)
#"""
ci.init_lens_sample( filename=lens_file, 
                     ntomo_bins=int(lens_ntomo))

ci.init_IA( ia_model = int(IA_model), 
            ia_redshift_evolution = int(IA_redshift_evolution))#"""

R, DS = DeltaSigma()

print(R, DS[0:30])
