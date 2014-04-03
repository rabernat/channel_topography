from pylab import *
import bump_channel
import os
import mycolors

# plot settings
rcParams['image.cmap'] = 'sun'
figdir = './figures/'

# the list of runs
runs = ['taux0125_rb0110_bump',
        'taux0250_rb0110_bump',
        'taux0500_rb0110_bump',
        'taux1000_rb0110_bump',
        'taux2000_rb0110_bump',
        'taux4000_rb0110_bump',
        'taux8000_rb0110_bump',
        'taux1000_rb0014_bump',
        'taux1000_rb0028_bump', 
        'taux1000_rb0055_bump',
        'taux1000_rb0220_bump',
        'taux1000_rb0440_bump',
        'taux1000_rb0880_bump',
        'taux0125_rb0110_flat',
        'taux0250_rb0110_flat',
        'taux0500_rb0110_flat',
        'taux1000_rb0110_flat',
        'taux2000_rb0110_flat',
        'taux4000_rb0110_flat',
        'taux8000_rb0110_flat',
        'taux2000_rb0110_bumplong'
         ]
N = len(runs)

# exclude the taux=8000 runs because they are still spinning up
#idx_wind = array([0,1,2,3,4,5,6])
#idx_flat = arange(13,N)
idx_wind = arange(7)
idx_flat = arange(13,N-1)

idx_rb = array([7,8,9,3,10,11,12])
# where to find the data

# paola's data
data_qg = {
  'Tau': array([0.0125,0.0250,0.0500,0.1000,0.2000,0.4000,0.8000]),
  'U2': array([0.0006,0.0011,0.0019,0.0034,0.0062,0.0118,0.0234]),
  'H1': array([0.4900,0.5477,0.6143,0.6988,0.8115,0.9693,1.2161]),
  'K':  array([0.1787,0.2777,0.4350,0.6846,1.0809,1.7100,2.7086]),
  'H1_nt': array([0.3887,0.5001,0.6386,0.8115,1.0279,1.2995,1.6409])
}

d = np.load('data/GAV_data.npz')
tau = d['GAV_tau0'][idx_wind]
r = d['GAV_rb'][idx_rb] / 2985. * 24*60*60 # express timescale in days

close('all')
#####################
# Big Global Figure #
#####################

rcParams['legend.fontsize'] = 7
rcParams['font.size'] = 7.
rcParams['lines.markersize'] = 4
fline = 'ok-'
rline = '^k--'
figure(figsize=(6.8,5.))

## stratification ##
#figure()
subplot(221)
plot(tau, d['GAV_D'][idx_flat], fline,
     tau, d['GAV_D'][idx_wind], rline,
#     data_qg['Tau'], data_qg['H1_nt']*1e3,'k-',
#     data_qg['Tau'], data_qg['H1']*1e3,'k--',
     )
plot(0.2, d['GAV_D'][-1], 'k*', markersize=10)
xlabel(r'$\tau_0$ (N m$^{-2}$)')
ylabel('h (m)')
yticks(range(800,1601,200))
legend(['flat', 'ridge'], loc='upper left')
title('Stratification Depth vs. Wind Stress')
grid()

# thermocline heat transport
subplot(222)
plot(tau, -d['GAV_VTD'][idx_flat], fline,
     tau, -d['GAV_VTD'][idx_wind], rline,)
plot(0.2, -d['GAV_VTD'][-1], 'k*', markersize=10)
xlabel(r'$\tau_0$ (N m$^{-2}$)')
ylabel(r"$-\widetilde{v_g \theta}$ (K m s$^{-1}$)")
title(r"$-\widetilde{v_g \theta}$ vs. Wind Stress")
grid()

# APE 
subplot(223)
plot(tau, d['GAV_APE'][idx_flat]/1e15, fline,
     tau, d['GAV_APE'][idx_wind]/1e15, rline)
plot(0.2, d['GAV_APE'][-1]/2/1e15, 'k*', markersize=10)
xlabel(r'$\tau_0$ (N m$^{-2}$)')
ylim([0,6])
ylabel(r'APE (PJ)')
grid()
title('Available Potential Energy vs. Wind Stress')
tight_layout()

# bottom vel
subplot(224)
loglog(tau, d['GAV_Ubotcen'][idx_flat], 'ok',
       tau, d['GAV_Ubotcen'][idx_wind], '^k',
       data_qg['Tau'], data_qg['U2'], 'k-')
loglog(0.2, d['GAV_Ubotcen'][-1], 'k*', markersize=10)
xlim([1e-2, 1e0])
xlabel(r'$\tau_0$ (N m$^{-2}$)')
ylabel(r'$\langle u \rangle_b$ (m/s)')
ylim([0.5e-3, 1])
grid()
title('Bottom Velocity vs. Wind Stress')

show()
savefig('%s/global_vars_APE.pdf' % figdir)

#################### END OF GLOBAL FIG ###################

