from pylab import *

HT = ma.masked_invalid(np.load('data/heat_transport_full_new.npz')['HT'])
HT_eqlat = ma.masked_invalid(np.load('data/heat_transport_full_new.npz')['HT_eqlat'])
tau = np.load('data/heat_transport_full.npz')['tau']
runtypes = ['flat', 'ridge']
pathnames = ['lat', r'$\Theta$', r'$\Psi$']
ssnames = ['','\Theta', '\Psi']
nicenames = ['Lat', 'Theta', 'Psi']
components = ['Mean','Standing','Transient']

close('all')
rcParams['font.size'] = 8
rcParams['legend.fontsize'] = 7
rcParams['grid.linewidth'] = 0.2
rcParams['lines.linewidth'] = 1.0
rcParams['lines.markersize'] = 4
for npath in arange(3):
    #figure(figsize=(6.8,2.5))
    figure(figsize=(3.5,5.))
    #for nrun in arange(2):
    for nrun in [1]:
        H_Ek_mean, H_SE_mean, H_TE_mean = zeros(len(tau)), zeros(len(tau)), zeros(len(tau))
        ratio, ratio_TE = zeros(len(tau)), zeros(len(tau))
        for n in arange(len(tau)):
            jr = (HT_eqlat[nrun,npath,n] > 1e6) & (HT_eqlat[nrun,npath,n] < 1.9e6)
            H_Ek = HT[nrun, npath, n, 0, :]
            H_SE = HT[nrun, npath, n, 1, :]
            H_TE = HT[nrun, npath, n, 2, :]
            H_Ek_mean[n] = H_Ek[jr].mean()
            H_SE_mean[n] = H_SE[jr].mean()
            H_TE_mean[n] = H_TE[jr].mean()
            ratio[n] = -(H_Ek/(H_SE+H_TE))[jr].mean()
            ratio_TE[n] = (H_TE/(H_SE+H_TE))[jr].mean()

        #subplot(1,2,mod(nrun+1,2))
        subplot(2,1,1)
        plot(tau, H_Ek_mean / 1e12, 'sk-',
             tau, H_SE_mean / 1e12, '^k-.',
             tau, H_TE_mean / 1e12, 'ok--')
        xlabel(r'$\tau_0$ (N m$^{-2}$)')
        ylabel(r'(TW)')
        legend([r'$\mathcal{H}_{Ek}^{%s}$' % ssnames[npath],
                r'$\mathcal{H}_{SE}^{%s}$' % ssnames[npath],
                r'$\mathcal{H}_{TE}^{%s}$' % ssnames[npath]], loc='upper left')
        title('Heat Transport vs. Wind: %s, %s contours' % (runtypes[nrun],pathnames[npath]))
        ylim([-400,400])
        grid()
        
        subplot(2,1,2)
        plot(tau, ratio, 'k-',
             tau, ratio_TE, 'k:')
        legend([r'$-\mathcal{H}_{Ek}^{%s}/(\mathcal{H}_{SE}^{%s}+\mathcal{H}_{TE}^{%s})$' %
                    (ssnames[npath],ssnames[npath],ssnames[npath]),
                r'$\mathcal{H}_{TE}^{%s}/(\mathcal{H}_{SE}^{%s}+\mathcal{H}_{TE}^{%s})$' %
                    (ssnames[npath],ssnames[npath],ssnames[npath])], loc='upper right')
        ylim([0,1.5])
        xlabel(r'$\tau_0$ (N m$^{-2}$)')
        title('Fractional Heat Transport: %s, %s contours' % (runtypes[nrun],pathnames[npath]))
        tight_layout()
        #savefig('%s/max_heat_transport-%s-%s.pdf' % (figdir, runtypes[nrun], pathnames[npath]))
        savefig('./figures/paper/mean_heat_transport-%s.pdf' % (nicenames[npath]))