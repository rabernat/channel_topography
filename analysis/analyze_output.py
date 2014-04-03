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
base_dir = os.path.join(os.environ['D'],'DATASTORE.RPA','projects','drag_strat')

# global average variables for each simulations
GAV_tau0 = empty(N)
GAV_rb = empty(N)
GAV_Utrans = empty(N)           # total zonal transport above 2000 m
GAV_Utrans_TW = empty(N)        # transport about 2000m due to thermal wind velocity
GAV_Ubarbot = empty(N)          # the mean zonally-averaged zonal velocity at z=-H
GAV_Ubotmax = empty(N)          # the max   "       "
GAV_Ubotcen = empty(N)          # average over middle of domain
GAV_Ubot = empty(N)             # the mean zonal velocity at the true bottom
GAV_Lmix = empty(N)             # mixing length
GAV_D = empty(N)                # thermocline depth
GAV_VTD = empty(N)              # vertically averaged heat flux within the thermocline
GAV_EKE_bot = empty(N)          # bottom EKE
GAV_EKE = empty(N)
GAV_APE = empty(N)              # available potential energy
GAV_Tp2 = empty(N)              # total transient eddy variance
GAV_work = empty((N,3))   # terms in the EKE budget
GAV_Uwork = empty((N,4))   # terms in the EKE budget
GAV_PE_budget = empty((N,2))   # terms in the EKE budget
GAV_VT_mean = empty(N)          # HT by mean circ
GAV_VT_standing = empty(N)      # HT by standing eddies
GAV_VT_transient = empty(N)     # HT by transient eddies
GAV_psi_iso = empty(N)          # residual circ
GAV_psi_bar = empty(N)          # mean circ
GAV_psi_eddy= empty(N)          # full eddy
GAV_psi_eddy_se = empty(N)      # standing
GAV_psi_eddy_te = empty(N)      # transient


#EKEbudget_legend = ['Useful Wind Work', 'Eddy Bottom Dissipation', 'Conversion', 'Net']
work_legend = [r'$\langle \tau \cdot u_g \rangle$',
               r"$r_b \langle |u_b - \bar u_b|^2 \rangle$",
               r"$\langle wb \rangle$", 'Net']
Uwork_legend = [r'$\langle \tau \cdot u_{TW} \rangle$',
                    r"$r_b \langle |u_b - \bar u_b|^2 \rangle$",
#                    r"$r_b \langle \bar v_b^2 \rangle$",
                    r"$\langle p_b H_x \bar u_b \rangle$",
                    r'$\langle wb \rangle$', 'Net']
PEbudget_legend = [r'$\langle wb \rangle$', r'$\langle \kappa b_z \rangle$']

Sv = 1e6

rcParams['figure.figsize'] = [6, 3.5]

#for n in [4,17]:
for n in arange(N):
    r = runs[n]
    try:
        b.reinit(r)
    except NameError:
        b = bump_channel.BumpChannel(base_dir, r)
    GAV_tau0[n] = b.tau0
    GAV_rb[n] = b.rb

    T = b.rdmds('Ttave')
    TT = b.rdmds('TTtave')
    U = b.rdmds('Uveltave')
    V = b.rdmds('Vveltave')
    Ubar = U.mean(axis=2)
    Vbar = V.mean(axis=2)
    Tbar = T.mean(axis=2)
    Ub = b.value_at_bottom(U)
    Vb = b.value_at_bottom(V)
    UU = b.rdmds('UUtave')
    VV = b.rdmds('VVtave')
    EKE = (UU-U*U) + (VV-V*V)
    UUb = b.value_at_bottom(UU)
    VVb = b.value_at_bottom(VV)

    GAV_EKE[n] = b.average_vertical(EKE).mean()
    GAV_Tp2[n] = b.average_vertical(TT - T**2).mean()
    GAV_APE[n] = b.get_APEcomputer().calc_APE()

    Vwave = V - Vbar[:,:,newaxis]
    Uwave = U - Ubar[:,:,newaxis]
    Twave = T - Tbar[:,:,newaxis] 
 
    # heat transport by zonal-averaged mean flow
    VbarTbar = b.integrate_vertical(Vbar*Tbar)
    # heat transport by the standing eddies
    #VTstanding = b.integrate_vertical((V*T).mean(axis=2) - Vbar*Tbar)
    VTstanding = b.integrate_vertical((Vwave*Twave).mean(axis=2))
    VTstanding_bt = (b.average_vertical(Vwave)*b.average_vertical(Twave)*b.depth).mean(axis=1)
    # heat transport by the fluctuating eddies
    VpTpfull = b.rdmds('VTtave').mean(axis=2) - (V*T).mean(axis=2)
    VpTpg = VpTpfull + (Vwave*Twave).mean(axis=2)
    VpTp = b.integrate_vertical(VpTpfull)
    VTnet = VbarTbar + VTstanding + VpTp
    VTlegend = ['zonal mean', 'standing', 'transient']
    
    
    # average values btw y=1000:1500 km
    jr = r_[0.5*b.Ny:0.75*b.Ny].astype('int')
    GAV_VT_transient[n] = VpTp[jr].mean()
    GAV_VT_standing[n] = VTstanding[jr].mean()
    GAV_VT_mean[n] = VbarTbar[jr].mean()
    
    try:
        psi_iso_z = b.get_psi_iso_z()
        psi_bar = b.get_psi_bar()
        lc = b.get_layers_computer()
        vflux,hv = lc.compute_vflux();
        psi_iso_se = -cumsum(vflux.mean(axis=2), axis=0) * b.Lx
        psi_iso_se_z = lc.transform_g_to_z(psi_iso_se,hv.mean(axis=2))
        psi_eddy_full = psi_iso_z - psi_bar
        psi_eddy_se = psi_iso_se_z - psi_bar
    except IOError:
        # some of them are missing
        psi_iso_z = zeros((b.Nz, b.Ny))
        psi_bar = psi_iso_z
        psi_eddy_full = psi_iso_z
        psi_eddy_se = psi_iso_z
    psi_eddy_te = psi_eddy_full - psi_eddy_se
    
    GAV_psi_bar[n] = psi_bar[16:].max()
    GAV_psi_eddy[n] = psi_eddy_full[16:].min()
    GAV_psi_eddy_se[n] = psi_eddy_se[16:].min()
    GAV_psi_eddy_te[n] = psi_eddy_te[16:].min()


    ## Zonal transport ##
    krange = r_[0:30]
    depth = sum(b.dzf[krange])
    
    # thermal wind velocities
    dbdy = b.ddy_cgrid_centered(b.g * b.tAlpha * T)
    u_tw = -cumsum( (dbdy*b.hFacW*b.dzf_grid)[::-1], axis=0)[::-1] / b.f[:,newaxis]
    # ssh velocity
    eta = b.rdmds('ETAtave')
    detady = b.ddy_cgrid_centered(eta)
    u_eta = -b.g * detady / b.f[:,newaxis]
    
    GAV_Utrans[n] = b.integrate_vertical(U,point='W',krange=krange).mean() * b.Ly
    GAV_Utrans_TW[n] = b.integrate_vertical(u_tw,point='W',krange=krange).mean() * b.Ly
    #GAV_Utrans_BT[n] = (Ub * b.depth).mean() * b.Ly

    GAV_Ubot[n] = b.value_at_bottom(U).mean()
    GAV_Ubarbot[n] = Ubar[-1].mean()
    GAV_Ubotcen[n] = Ubar[-1,100:300].mean()
    GAV_Ubotmax[n] = Ubar[-1].max()

    # thermocline depth
    # only integrate down to top of ridge
    D = -2 * b.integrate_vertical( b.zc[:,newaxis,newaxis]
        * (T - b.value_at_bottom(T)), krange=krange) /  b.integrate_vertical(T,krange=krange)
    GAV_D[n] = D[-1].mean()
    
    # average over top D meters and the northern half of the domain
    Dkrange = find(b.zc > -GAV_D[n])
    GAV_VTD[n] = b.average_vertical(VpTpg, krange=Dkrange)[b.Ny/2:].mean()
        
    # mixing length
    Tp2 = b.rdmds('TTtave') - T**2
    dTdy_bulk = 8. / (2000e3)
    GAV_Lmix[n] = Tp2[:15].mean() / dTdy_bulk
    
    ## TERMS OF ZONAL/VERTICAL MOMENTUM BUDGET ##

    # 0 - wind stress
    # 1 - form drag
    # 2 - bottom friction
    # 3 - advective flux divergence
    # the equations
    mom_bud_legend = [
        'wind stress', 'form stress', 'bottom friction', 'advection', 'total'
    ]
    mom_bud_values = empty((5,b.Ny))

    # b.wind_stress already exists
    mom_bud_values[0] = b.wind_stress

    # form drag
    pbot = b.rdmds('PHLtave')
    # centered difference
    dHdx = -0.5 * (b.depth[:,r_[1:b.Nx,0]] - b.depth[:,r_[b.Nx-1,:b.Nx-1]]) / b.dxc[0]
    form_stress = b.rho0 * (pbot*dHdx).mean(axis=1)
    mom_bud_values[1] = -form_stress

    # frictional bottom drag
    mom_bud_values[2] = -b.rho0 * b.rb * Ub.mean(axis=1)

    # advective momentum flux divergence (includes eddy and mean)
    UV = b.integrate_vertical(b.rdmds('UVtave')).mean(axis=1)
    # divergence (defined at V point)
    dUVdy = (UV - UV[r_[1:b.Ny,b.Ny-1]]) / b.dyc[0]
    mom_bud_values[3] = b.rho0 * dUVdy

    # sum
    mom_bud_values[4] = sum(mom_bud_values[:4],axis=0)

    N2 = b.g * b.tAlpha * b.ddz_cgrid_centered(Tbar)

    ## TERMS OF ENERGY BUDGET ##

    # 0 - wind work on geostrophic surface flow
    # 1 - bottom drag dissipation
    # 2 - conversion
    
    # units are W m^2
    # this is the same as doing a global integral (3D) of the energy budget
    # and dividing by the surface area
    
    # calculate the surface geostrophic U using the mean ssh
    eta = b.rdmds('ETAtave')
    detady = b.ddy_cgrid_centered(eta)
    Us_g = -b.g * detady / b.f[:,newaxis]

    # calculate the surface geostrophic U by integrating thermal wind
    Us_tw = u_tw[0]
    
    # conversion term
    WB = b.g * b.tAlpha * b.integrate_vertical(b.rdmds('WTtave'))
    Bdif = b.g * b.tAlpha * b.integrate_vertical(b.rdmds('Tdiftave'))
    
    GAV_work[n,0] = mean(b.wind_stress * Us_g.mean(axis=1))
    GAV_work[n,1] = -b.rho0 * b.rb * (UUb + VVb).mean()
    GAV_work[n,2] = b.rho0 * WB.mean()
    
    GAV_Uwork[n,0] = mean(b.wind_stress * Us_tw.mean(axis=1)) 
    GAV_Uwork[n,1] = -b.rho0 * b.rb * (UUb - Ub.mean(axis=1)[:,newaxis] + VVb).mean()
    GAV_Uwork[n,2] = - (form_stress * Ub.mean(axis=1)).mean()
    GAV_Uwork[n,3] = b.rho0 * WB.mean()
    
    GAV_PE_budget[n,0] = b.rho0 * WB.mean()
    GAV_PE_budget[n,1] = b.rho0 * Bdif.mean()
    
    ## PLOTS ##
    close('all')
    tit = r'$\tau_0 =$ %4.3f N m$^{-2}$, $r =$ %4.3f days$^{-1}$ (%s)' % (b.tau0, b.rb / 2985. * 24*60*60, b.bathy)
    # cheat to make y prettier
    Y = linspace(0, b.Ly/1000., b.Ny)
    
    if False:
        figure()
        plot(Y, mom_bud_values[:4].T)
        plot(Y, mom_bud_values[4], 'k', linewidth=2)
        grid()
        xlabel('Y (km)')
        ylabel('Stress (N m$^{-2}$)')
        title('Zonal Momentum Budget : %s' % tit)
        legend(mom_bud_legend)
        tight_layout()
        savefig('%s/runs/momentum_budget-%s.pdf' % (figdir, r))
    
        figure()
        Qfac = b.rho0 * 3994 * b.Lx / 1e12
        plot(Y, Qfac*VbarTbar, 'k', Y, Qfac*VTstanding, 'm',
            Y, Qfac*VpTp, 'c', Y, Qfac*VTnet, '0.7')
        grid()
        xlabel('Y (km)')
        ylabel('Heat Transport (TW)')
        if b.tau0==0.2:
            ylim([-120,120])
        legend(VTlegend, 'upper left')
        title(tit)
        savefig('%s/runs/heat_transport-%s.pdf' % (figdir, r))
    
    # zonal average plots
    if False:
        Ulevs = arange(0,13)
        EKElevs = arange(0,51,5)*1e-3
        Tlevs = arange(0,7.1)
        psi_levs = arange(-10,11)/10.
    
        figure()
        contourf(Y, b.zc, psi_iso_z / Sv, psi_levs, cmap=get_cmap('posneg'), extend='both')
        xlabel('Y (km)'); ylabel('Z (m)')
        title('$\Psi_{iso}$ (Sv) : ' + tit)
        tight_layout()
        colorbar()
        contour(Y, b.zc, Tbar, Tlevs, colors='k')
        savefig('%s/runs/psi_iso-%s.pdf' % (figdir, r))

        set_cmap(get_cmap('posneg'))
        figure()
        psi_levs = arange(-6,6)+0.5
        contourf(Y,b.zc,psi_bar / Sv, psi_levs, extend='both')
        xlabel('Y (km)'); ylabel('Z (m)')
        title(r'$\bar \Psi$ (Sv) : ' + tit)
        tight_layout()
        colorbar(ticks=arange(-5,6))
        contour(Y, b.zc, Tbar, Tlevs, colors='k')
        savefig('%s/runs/psi_bar-%s.pdf' % (figdir, r))

        figure()
        psi_levs = arange(-6,6)+0.5
        contourf(Y,b.zc,psi_eddy_full /Sv, psi_levs, extend='both')
        xlabel('Y (km)'); ylabel('Z (m)')
        title(r'$\Psi^\ast$ (Sv) : ' + tit)
        tight_layout()
        colorbar(ticks=arange(-5,6))
        contour(Y, b.zc, Tbar, Tlevs, colors='k')
        savefig('%s/runs/psi_eddy_full-%s.pdf' % (figdir, r))

        figure()
        psi_levs = arange(-6,6)+0.5
        contourf(Y,b.zc,psi_eddy_se / Sv, psi_levs, extend='both')
        xlabel('Y (km)'); ylabel('Z (m)')
        title(r'$\Psi^\ast_{SE}$ (Sv) : ' + tit)
        tight_layout()
        colorbar(ticks=arange(-5,6))
        contour(Y, b.zc, Tbar, Tlevs, colors='k')
        savefig('%s/runs/psi_eddy_se-%s.pdf' % (figdir, r))

        figure()
        psi_levs = arange(-6,6)+0.5
        contourf(Y,b.zc,psi_eddy_te / Sv, psi_levs, extend='both')
        xlabel('Y (km)'); ylabel('Z (m)')
        title(r'$\Psi^\ast_{TE}$ (Sv) : ' + tit)
        tight_layout()
        colorbar(ticks=arange(-5,6))
        contour(Y, b.zc, Tbar, Tlevs, colors='k')
        savefig('%s/runs/psi_eddy_te-%s.pdf' % (figdir, r))
        # figure()
        # pcolormesh(Y, b.zc, Kb, cmap=get_cmap('jet'))
        # xlabel('Y (km)'); ylabel('Z (m)')
        # title('$K_b$ (m$^2$s$^{-1}$) : ' + tit)
        # tight_layout()
        # clim([0,8000])
        # colorbar()
        # contour(Y, b.zc, Tbar, Tlevs, colors='k')
        # savefig('%s/runs/Kb-%s.pdf' % (figdir, r))
    
        figure()
        PVlevs = arange(0,15)*1e-10
        contourf(Y, b.zc, -b.f0 * N2, PVlevs, extend='max')
        xlabel('Y (km)'); ylabel('Z (m)')
        title('PV (Sv) : ' + tit)
        tight_layout()
        colorbar()
        contour(Y, b.zc, Tbar, Tlevs, colors='k')
        savefig('%s/runs/PV-%s.pdf' % (figdir, r))
 
        figure()
        # contourf(Y, b.zc, EKE.mean(axis=2), EKElevs, extend='max')
        c=contourf(Y, b.zc, u_tw.mean(axis=2)*100, Ulevs, linewidth=2)
        xlabel('Y (km)'); ylabel('Z (m)')
        title(r'$\bar u_{tw}$ (ms$^{-1}$) : ' + tit)
        tight_layout()
        colorbar()
        contour(Y, b.zc, Tbar, Tlevs, colors='k')
        # text(300, -2000, 'zonal transport = %4.1f Sv' % (GAV_Utrans[n] / Sv) )
        # text(300, -2200, 'barotropic transport = %4.1f Sv' % (GAV_Utrans_BT[n] / Sv) )
        if b.bathy=='bump':
            gca().fill_between(Y[r_[0,-1]], -array([2058.0,2058.0]),-array([b.H,b.H]),
                    edgecolor='none', facecolor='0.7', alpha=0.3)
        savefig('%s/runs/Ubar-%s.pdf' % (figdir, r))
 
    
        # figure()
        # contourf(Y, b.zc, EKE.mean(axis=2), EKElevs, extend='max')
        # xlabel('Y (km)'); ylabel('Z (m)')
        # title('EKE (m$^2$s$^{-2}$) : ' + tit)
        # tight_layout()
        # colorbar()
        # c=contour(Y, b.zc, Ubar*100, Ulevs, colors='g', linewidth=2)
        # clabel(c,array([0.]), fmt='U = %2d')
        # contour(Y, b.zc, Tbar, Tlevs, colors='k')
        # text(300, -2000, 'zonal transport = %4.1f Sv' % (GAV_Utrans[n] / Sv) )
        # text(300, -2200, 'barotropic transport = %4.1f Sv' % (GAV_Utrans_BT[n] / Sv) )
        # savefig('%s/runs/EKE_U_T_section-%s.pdf' % (figdir, r))

GAV_Utrans_BT = GAV_Utrans - GAV_Utrans_TW

# save data
mydat = dict()
for k in locals().keys():
    if k[:3] == 'GAV':
        mydat[k] = locals()[k]
        print(k)
savez('data/GAV_data.npz', **mydat)


# paola's data
data_qg = {
  'Tau': array([0.0125,0.0250,0.0500,0.1000,0.2000,0.4000,0.8000]),
  'U2': array([0.0006,0.0011,0.0019,0.0034,0.0062,0.0118,0.0234]),
  'H1': array([0.4900,0.5477,0.6143,0.6988,0.8115,0.9693,1.2161]),
  'K':  array([0.1787,0.2777,0.4350,0.6846,1.0809,1.7100,2.7086]),
  'H1_nt': array([0.3887,0.5001,0.6386,0.8115,1.0279,1.2995,1.6409])
}


# global plots
close('all')
#rcParams['figure.figsize'] = [6., 4.]

tau = GAV_tau0[idx_wind]
r = GAV_rb[idx_rb] / 2985. * 24*60*60 # express timescale in days

for bathy in ['bump', 'flat']:
    if bathy=='bump':
        idx = idx_wind
    else:
        idx = idx_flat        

    ## transport ##
    Utrans_legend = ['Total', 'Baroclinic', 'Barotropic']
    figure()
    plot(tau, GAV_Utrans[idx]/Sv, 'k.-',
            tau, GAV_Utrans_TW[idx]/Sv, 'b.-',
            tau, GAV_Utrans_BT[idx]/Sv, 'r.-')
    #ylim([0,120])
    xlabel(r'$\tau_0$ (N m$^{-2}$)')
    ylabel('Transport (Sv)')
    legend(Utrans_legend, loc='upper left')
    title('Zonal Transport vs. Wind Stress (%s)' % bathy)
    grid()
    tight_layout()
    savefig('%s/transport_vs_wind_%s.pdf' % (figdir, bathy))

    ## Energy Budget ##
    figure()
    plot(tau, GAV_work[idx], '.-')
    plot(tau, GAV_work[idx].sum(axis=1), 'k.-')
    legend(work_legend, loc='upper left')
    xlabel(r'$\tau_0$ (N m$^{-2}$)')
    ylabel(r'Conversion Rate (W m$^{-2}$)')
    title('Kinetic Energy Sources / Sinks (%s)' % bathy)
    #ylim([-0.02,0.02])
    grid()
    tight_layout()
    savefig('%s/work_%s.pdf' % (figdir, bathy))
    
    ## Useful Energy Budget ##
    figure()
    plot(tau, GAV_Uwork[idx], '.-')
    plot(tau, GAV_Uwork[idx].sum(axis=1), 'k.-')
    legend(Uwork_legend, loc='upper left')
    xlabel(r'$\tau_0$ (N m$^{-2}$)')
    ylabel(r'Conversion Rate (W m$^{-2}$)')
    title('Useful Kinetic Energy Sources / Sinks (%s)' % bathy)
    #ylim([-0.02,0.02])
    grid()
    tight_layout()
    savefig('%s/useful_work_%s.pdf' % (figdir, bathy)) 
       
    figure()
    plot(tau, GAV_PE_budget[idx], '.-')
    plot(tau, GAV_PE_budget[idx].sum(axis=1), 'k.-')
    legend(PEbudget_legend, loc='upper left')
    xlabel(r'$\tau_0$ (N m$^{-2}$)')
    ylabel(r'Conversion Rate (W m$^{-2}$)')
    title('Potential Energy Sources / Sinks (%s)' % bathy)
    grid()
    tight_layout()
    savefig('%s/PE_budget_wind_%s.pdf' % (figdir, bathy))

# bulk diffusivity
dtdy_mean = 8 / 2e6
Knet0 = GAV_tau0 * 2 / pi * b.Ly / (b.rho0 * abs(b.f0) * GAV_D)
Knet = - (GAV_VT_transient + GAV_VT_standing) / (GAV_D*dtdy_mean)
Kstand = - GAV_VT_standing / (GAV_D*dtdy_mean)
Ktrans = - GAV_VT_transient / (GAV_D*dtdy_mean)

# log-log diff
figure()
loglog(tau, Knet[idx_flat] / GAV_EKE[idx_flat]**0.5, 'ok-',
      tau, Knet[idx_wind] / GAV_EKE[idx_wind]**0.5, 'oc-',)
xlabel(r'$\tau_0$ (N m$^{-2}$)')
ylabel(r'$L_{eff}$ (m)')
legend(['flat', 'ridge'], loc='upper left')
title('Bulk Mixing Length vs. Wind Stress')
grid()
tight_layout()
savefig('%s/Lbulk_vs_wind.pdf' % figdir)

#####################
# Big Global Figure #
#####################

#rc('figure.subplot', left=0.2, right=0.95, bottom=0.1, top=0.95, wspace=0.1, hspace=0.25)
rcParams['legend.fontsize'] = 7
rcParams['font.size'] = 7.
rcParams['lines.markersize'] = 4
fline = 'ok-'
rline = '^k--'
# rline ='oc-'
# diff
figure(figsize=(6.8,5.))
#figure(figsize=(3.3,6.))
#subplot(223)
#plot(tau, Knet[idx_flat], fline,
#     tau, Knet[idx_wind], rline,)
#gca().fill_between(tau, Knet[0]* (tau/tau[0]), 
#                     Knet[idx_flat][0] * (tau/tau[0])**0.5,
#                     edgecolor='none', color='0.8')
#xlabel(r'$\tau_0$ (N m$^{-2}$)')
#ylabel(r'$K_g$ (m$^2$ s$^{-1}$)')
#legend(['flat', 'ridge'], loc='upper left')
#title('$K_g$ vs. Wind Stress')
#grid()
#ylim([0,10000])
#tight_layout()
#savefig('%s/Kg_vs_wind.pdf' % figdir)


## stratification ##
#figure()
subplot(221)
plot(tau, GAV_D[idx_flat], fline,
     tau, GAV_D[idx_wind], rline,
#     data_qg['Tau'], data_qg['H1_nt']*1e3,'k-',
#     data_qg['Tau'], data_qg['H1']*1e3,'k--',
     )
plot(0.2, GAV_D[-1], 'k*', markersize=10)
xlabel(r'$\tau_0$ (N m$^{-2}$)')
ylabel('h (m)')
yticks(range(800,1601,200))
legend(['flat', 'ridge'], loc='upper left')
title('Stratification Depth vs. Wind Stress')
grid()
#tight_layout()
#savefig('%s/stratification_vs_wind.pdf' % figdir)

# log version
# figure()
# loglog(tau, GAV_D[idx_flat], 'ok-',
#      tau, GAV_D[idx_wind], 'oc-',)
# xlabel(r'$\tau_0$ (N m$^{-2}$)')
# ylabel('h (m)')
# legend(['flat', 'ridge'], loc='upper left')
# title('Stratification Depth vs. Wind Stress')
# grid()
# tight_layout()
# savefig('%s/stratification_vs_wind_loglog.pdf' % figdir)

## thermal wind transport ##
#figure()
#subplot(222)
#plot(tau, GAV_Utrans_TW[idx_flat]/1e6, fline,
#     tau, GAV_Utrans_TW[idx_wind]/1e6, rline,)
#xlabel(r'$\tau_0$ (N m$^{-2}$)')
#ylabel('T (Sv)')
#legend(['flat', 'ridge'], loc='upper left')
#title('Thermal Wind Transport vs. Wind Stress')
#grid()
#tight_layout()
#savefig('%s/tw_transport_vs_wind.pdf' % figdir)

# thermocline heat transport
subplot(222)
plot(tau, -GAV_VTD[idx_flat], fline,
     tau, -GAV_VTD[idx_wind], rline,)
plot(0.2, -GAV_VTD[-1], 'k*', markersize=10)
xlabel(r'$\tau_0$ (N m$^{-2}$)')
ylabel(r"$-\widetilde{v_g \theta}$ (K m s$^{-1}$)")
#legend(['flat', 'ridge'], loc='upper left')
title(r"$-\widetilde{v_g \theta}$ vs. Wind Stress")
grid()

#layer = np.load('data/U2_vs_tau_2layer.npz')
# bottom vel
subplot(224)
loglog(tau, GAV_Ubotcen[idx_flat], 'ok',
       tau, GAV_Ubotcen[idx_wind], '^k',
       data_qg['Tau'], data_qg['U2'], 'k-')
loglog(0.2, GAV_Ubotcen[-1], 'k*', markersize=10)
#loglog([1e-3,1e1],[0.2e-3, 0.2e1],':k')
xlim([1e-2, 1e0])
xlabel(r'$\tau_0$ (N m$^{-2}$)')
ylabel(r'$\bar u$ (m/s)')
ylim([0.5e-3, 1])
grid()
#legend(['ridge (mean)', 'flat (mean)', 'ridge (center)', 'flat (center)'], loc='lower right')
title('Bottom Velocity vs. Wind Stress')
#savefig('figures/Ubot_vs_wind_loglog.pdf')

# EKE #
#figure()
subplot(223)
loglog(tau, GAV_EKE[idx_flat], fline,
       tau, GAV_EKE[idx_wind], rline)
gca().fill_between(tau, GAV_EKE[idx_flat][0]*1.2 * (tau/tau[0]), 
                        GAV_EKE[idx_flat][0]*0.95 * (tau/tau[0])**0.5,
                        edgecolor='none', color='0.8')
loglog(0.2, GAV_EKE[-1], 'k*', markersize=10)
ylim([1e-3, 1.4e-1])
xlabel(r'$\tau_0$ (N m$^{-2}$)')
ylabel(r'EKE (m$^2$ s$^{-2}$)')
#legend(['flat', 'ridge'], loc='upper left')
grid()
title('Eddy Kinetic Energy vs. Wind Stress')
tight_layout()
#savefig('%s/EKE_vs_wind_loglog.pdf' % figdir)
#savefig('%s/global_vars_vs_wind.pdf' % figdir)
savefig('%s/global_vars_latest.pdf' % figdir)

#################### END OF GLOBAL FIG ###################

figure()
plot(r, GAV_Utrans[idx_rb]/Sv, 'k.-',
        r, GAV_Utrans_TW[idx_rb]/Sv, 'b.-',
        r, GAV_Utrans_BT[idx_rb]/Sv, 'r.-')
ylim([0,120])
xlabel('$r$ (days$^{-1}$)')
ylabel('Transport (Sv)')
grid()
#legend(Utrans_legend, loc='upper left')
title('Zonal Transport vs. Bottom Drag')
tight_layout()
savefig('%s/transport_vs_drag.pdf' % figdir)


#figure()
#plot(r, GAV_EKE[idx_rb], 'k.-',
#        r, GAV_EKE_bot[idx_rb], 'r.-')
#ylim([0,0.01])
#xlabel('$r$ (days$^{-1}$)')
#ylabel('EKE (m$^2$s$^{-2}$)')
#legend(EKE_legend, loc='upper left')
#title('EKE vs. Bottom Drag')
#grid()
#tight_layout()
#savefig('%s/EKE_vs_drag.pdf' % figdir)

#figure()
#plot(r, GAV_EKE_budget[idx_rb,:3], '.-')
#plot(r, GAV_EKE_budget[idx_rb,3], 'k.-')
#legend(EKEbudget_legend, loc='upper left')
#xlabel('$r$ (days$^{-1}$)')
#title('Kinetic Energy Sources / Sinks')
#grid()
#tight_layout()
#savefig('%s/KE_budget_drag.pdf' % figdir)

figure()
plot(r, GAV_PE_budget[idx_rb,:2], '.-')
plot(r, GAV_PE_budget[idx_rb,2], 'k.-')
legend(PEbudget_legend, loc='upper left')
xlabel('$r$ (days$^{-1}$)')
title('Potential Energy Sources / Sinks')
grid()
tight_layout()
savefig('%s/PE_budget_drag.pdf' % figdir)


# heat transport figure
figure()
plot(tau, GAV_VT_mean[idx_flat], '.k-')
plot(tau, GAV_VT_transient[idx_flat], '.c-')
xlabel(r'$\tau_0$ (N m$^{-2}$)')
legend(['Mean (flat)','Transient (flat)'], loc='upper left')
ylabel(r'Heat Transport (TW)')
grid()
tight_layout()
savefig('%s/HT_flat.pdf' % figdir)
plot(tau, GAV_VT_mean[idx_wind], '.k--')
plot(tau, GAV_VT_transient[idx_wind], '.c--')
plot(tau, GAV_VT_standing[idx_wind], '.m--')
legend(['Mean (flat)','Transient (flat)','Mean (ridge)','Transient (ridge)','Standing (ridge)'], loc='upper left')
savefig('%s/HT_all.pdf' % figdir)

# MOC figure
figure()
plot(tau, GAV_psi_bar[idx_flat] / Sv, '.k-')
plot(tau, GAV_psi_eddy[idx_flat] / Sv, '.c-')
xlabel(r'$\tau_0$ (N m$^{-2}$)')
legend([r'$\bar \Psi$ (flat)',r'$\Psi^\ast$ (flat)'], loc='upper left')
ylabel(r'Transport (Sv)')
grid()
tight_layout()
savefig('%s/psi_vs_tau_flat.pdf' % figdir)
plot(tau, GAV_psi_bar[idx_wind] / Sv, '.k--')
plot(tau, GAV_psi_eddy_te[idx_wind] /Sv, '.c--')
plot(tau, GAV_psi_eddy_se[idx_wind]/Sv, '.m--')
legend([r'$\bar \Psi$ (flat)',r'$\Psi^\ast$ (flat)',r'$\bar \Psi$ (ridge)',r'$\Psi^\ast_{TE}$ (ridge)', r'$\Psi^\ast_{SE}$ (ridge)'], loc='upper left')
savefig('%s/psi_vs_tau_all.pdf' % figdir)

# new HT figures
HT = ma.masked_invalid(np.load('data/heat_transport_full.npz')['HT'])
HT_eqlat = ma.masked_invalid(np.load('data/heat_transport_full.npz')['HT_eqlat'])
tau = np.load('data/heat_transport_full.npz')['tau']
# HT = zeros((len(runs), len(pathtypes), len(components), len(mytau)))
runtypes = ['flat', 'ridge']
pathnames = ['lat', r'$\Theta$', r'$\Psi$']
ssnames = ['','\Theta', '\Psi']
nicenames = ['Lat', 'Theta', 'Psi']
components = ['Mean','Standing','Transient']


rcParams['font.size'] = 8
rcParams['legend.fontsize'] = 7
rcParams['grid.linewidth'] = 0.2
rcParams['lines.linewidth'] = 1.0
rcParams['lines.markersize'] = 4
for npath in arange(3):
    figure(figsize=(6.8,2.5))
    for nrun in arange(2):
        #j = argmax(HT[nrun,npath,0],axis=1)
        j = argmin((HT_eqlat[nrun,npath]-1250e3)**2, axis=1)
        #plot(tau, HT[nrun, npath, 0, :], '.k-',
        #     tau, HT[nrun, npath, 1, :], '.m-',
        #     tau, HT[nrun, npath, 2, :], '.c-')
        subplot(1,2,mod(nrun+1,2))
        plot(tau, HT[nrun, npath, 0, arange(len(tau)), j], 'sk-',
             tau, HT[nrun, npath, 1, arange(len(tau)), j], '^k-.',
             tau, HT[nrun, npath, 2, arange(len(tau)), j], 'ok--')
        xlabel(r'$\tau_0$ (N m$^{-2}$)')
        ylabel(r'Max Heat Transport (TW)')
        legend([r'$\mathcal{H}_{EK}^{%s}$' % ssnames[npath],
                r'$\mathcal{H}_{SE}^{%s}$' % ssnames[npath],
                r'$\mathcal{H}_{TE}^{%s}$' % ssnames[npath]], loc='upper left')
        title('Heat Transport vs. Wind: %s, %s contours' % (runtypes[nrun],pathnames[npath]))
        ylim([-400,400])
        grid()
        tight_layout()
        #savefig('%s/max_heat_transport-%s-%s.pdf' % (figdir, runtypes[nrun], pathnames[npath]))
        savefig('./figures/paper/max_heat_transport-%s.pdf' % (nicenames[npath]))
        
figure(figsize=(3.2,2.8))
loglog(tau, GAV_Ubarbot[idx_wind], '^-k')
loglog(tau, GAV_Ubarbot[idx_flat], 'o-k')
loglog(tau, GAV_Ubotcen[idx_wind], '^--k')
loglog(tau, GAV_Ubotcen[idx_flat], 'o--k')
loglog([1e-3,1e1],[0.2e-4, 0.2e0],':k')
xlim([1e-2, 1e0])
xlabel(r'$\tau_0$ (N m$^{-2}$)')
ylabel(r'$\bar u$ (m/s)')
legend(['ridge (mean)', 'flat (mean)', 'ridge (center)', 'flat (center)'], loc='lower right')
title('Bottom Velocity vs. Wind')
tight_layout()
savefig('figures/Ubot_vs_wind_loglog.pdf')