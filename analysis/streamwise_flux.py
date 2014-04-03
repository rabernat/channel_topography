from pylab import *
import os
import bump_channel
import MITgcmdata.poisson
import mycolors

rcParams['figure.figsize'] = [6, 3.5]
rcParams['legend.fontsize'] = 7
rcParams['font.size'] = 8

base_dir = os.path.join(os.environ['D'],'projects','drag_strat')
figdir = './figures/runs_HT'

mytau = array([0.0125, 0.0250, 0.05, 0.1, 0.2, 0.4, 0.8])
runs = ['flat', 'bump']
pathtypes = ['LAT', 'PSI', 'THETA']
pathnames = ['Lat', r'$\Psi$', r'$\Theta$']
components = ['Mean','Standing','Transient']

# spacing in streamline space
Nt = 100

# heat transport arrays
HT = zeros((len(runs), len(pathtypes), len(mytau), len(components), Nt))
HT_coords = zeros((len(runs), len(pathtypes), len(mytau), Nt))
HT_eqlat = zeros((len(runs), len(pathtypes), len(mytau), Nt))

doVertical = True

#for nrun in arange(len(runs)):
#for nrun in [0,]:
for nrun in [1,]:
    run = runs[nrun]
#    for ntau in arange(len(mytau)):
    for ntau in [4,]:

        close('all')

        tau = mytau[ntau]
        r = 'taux%04d_rb0110_%s' % (tau*10000, run)
        b = bump_channel.BumpChannel(base_dir, r)
        sol = MITgcmdata.poisson.solver(b) # conjugate gradient solver

        T = b.rdmds('Ttave')
        Tw = b.cgrid_to_ugrid(T) # on the U grid
        Ts = b.cgrid_to_vgrid(T) # on the V grid
        V = b.rdmds('Vveltave')
        U = b.rdmds('Uveltave')
        P = b.rdmds('Phhytave')

        Ug = -(b.f[newaxis,:,newaxis]**-1) * b.ddy_cgrid_to_vgrid(P)
        Vg = (b.f[newaxis,:,newaxis]**-1) * b.ddx_cgrid_to_ugrid(P)
        Ua = U - Ug
        Va = V - Vg

        UT = b.rdmds('UTtave')
        VT = b.rdmds('VTtave')

        # vertically integrated fluxes

        THETA = b.average_vertical(T)
        THETAw = b.cgrid_to_ugrid(THETA)
        THETAs = b.cgrid_to_ugrid(THETA)
        
        # full time-mean flux
        divTflux = b.ddx_ugrid_to_cgrid(U*Tw) + b.ddy_vgrid_to_cgrid(V*Ts)

        # ageostrophic flux
        divTflux_a = b.ddx_ugrid_to_cgrid(Ua*Tw) + b.ddy_vgrid_to_cgrid(Va*Ts)

        # geostrophic flux
        divTflux_g = b.ddx_ugrid_to_cgrid(Ug*Tw) + b.ddy_vgrid_to_cgrid(Vg*Ts)

        # eddy flux
        divTflux_e = b.ddx_ugrid_to_cgrid(UT - U*Tw) + b.ddy_vgrid_to_cgrid(VT - V*Ts)

        #for npath in arange(len(pathtypes)):
        for npath in [2,]:
        
            avg = pathtypes[npath]
            
            # vertical average
            if avg=='LAT':
                T0 = tile(b.yc[:,newaxis], [1, b.Nx]) / 1e3
                t = linspace(0, b.Ly, Nt) / 1e3
            elif avg=='THETA':
                T0 = b.average_vertical(T)
                tmean = T0.mean(axis=1)
                t = linspace(tmean.min(), tmean.max(), Nt)
            elif avg=='PSI':
                T0 = cumsum(b.integrate_vertical(U), axis=0) * b.dyc[0] / 1e6
                tmean = T0.mean(axis=1)
                extra=0.
                if run=='bump':
                    extra = tmean.max() * log2(tau/mytau[0])/7.5
                t = linspace(tmean.min()-extra, tmean.max()+extra, Nt)    
   
            # output arrays
            if doVertical:
                Tflux_a_t = zeros((b.Nz,Nt))
                Tflux_g_t = zeros((b.Nz,Nt))
                Tflux_e_t = zeros((b.Nz,Nt))
            else:
                Tflux_a_t = zeros(Nt)
                Tflux_g_t = zeros(Nt)
                Tflux_e_t = zeros(Nt)
                divTflux_a_vint = b.integrate_vertical(divTflux_a)
                divTflux_g_vint = b.integrate_vertical(divTflux_g)
                divTflux_e_vint = b.integrate_vertical(divTflux_e)
            
            grad_T0 = (b.ddx_cgrid_centered(T0)**2 + b.ddy_cgrid_centered(T0)**2)**0.5
            area = zeros(Nt)
            aint_grad_T0 = zeros(Nt)
            aint_T_gradT0 = zeros((b.Nz,Nt))
            

            # divergence method
            tbnds = hstack([-inf, 0.5*(t[1:]+t[:-1]), inf])
            for n in arange(Nt):
                idx = T0<=tbnds[n+1]
                area[n] = b.rac[idx].sum()
                aint_grad_T0[n] = (b.rac*grad_T0)[idx].sum()
                if doVertical:
                    for k in arange(b.Nz):
                        aint_T_gradT0[k,n] = (b.rac*T[k]*grad_T0)[idx].sum()
                        Tflux_a_t[k,n] = (b.rac*divTflux_a[k])[idx].sum() 
                        Tflux_g_t[k,n] = (b.rac*divTflux_g[k])[idx].sum()
                        Tflux_e_t[k,n] = (b.rac*divTflux_e[k])[idx].sum()
                else:
                    Tflux_a_t[n] = (b.rac*divTflux_a_vint)[idx].sum() 
                    Tflux_g_t[n] = (b.rac*divTflux_g_vint)[idx].sum()
                    Tflux_e_t[n] = (b.rac*divTflux_e_vint)[idx].sum()

            Y0 = area / b.Lx
            Y00 = 0.5*(Y0[1:] + Y0[:-1])
            #Y0 = interp(t,T0.mean(axis=1), b.yc)
            # I can't believe this worked...
            if doVertical:
                Tbar_T0 = -diff(aint_T_gradT0) / (t[1]-t[2]) / (-diff(aint_grad_T0)/(t[1]-t[2]))[newaxis,:]

            # vertical integrals
            Qfac = b.rho0 * 3994
            if doVertical:
                H_Ek = Qfac*b.integrate_vertical(Tflux_a_t)
                H_SE = Qfac*b.integrate_vertical(Tflux_g_t)
                H_TE = Qfac*b.integrate_vertical(Tflux_e_t)
            else:
                H_Ek = Qfac*Tflux_a_t
                H_SE = Qfac*Tflux_g_t
                H_TE = Qfac*Tflux_e_t

            # HT = zeros((len(runs), len(pathtypes), len(components), len(mytau)))
            HT[nrun, npath, ntau, 0] = H_Ek
            HT[nrun, npath, ntau, 1] = H_SE
            HT[nrun, npath, ntau, 2] = H_TE
            HT_coords[nrun, npath, ntau] = t
            HT_eqlat[nrun, npath, ntau] = Y0
            
            # figures
            if doVertical:
                myrun = run
                if myrun=='bump':
                    myrun='ridge'
                if avg=='LAT':
                    xlab = 'y (km)'
                    hname = r'\mathcal{H}'
                    sname = 'y'
                    mytit = 'Latitude Circles - %s' % myrun
                elif avg=='THETA':
                    xlab = r'$\Theta$ ($^\circ$ C)'
                    hname = r'\mathcal{H}^\Theta'
                    sname = '\Theta'
                    mytit = r'$\Theta$ Contours - %s' % myrun
                elif avg=='PSI':
                    xlab = r'$\Psi$ (Sv)'        
                    hname = r'\mathcal{H}^\Psi'
                    sname = '\Psi'
                    mytit = r'$\Psi$ Contours - %s' % myrun
            
                # close('all')
                # rcParams['font.size'] = 8.
                # rc('figure.subplot', left=0.2, right=0.93, bottom=0.13, top=0.94,
                #                     hspace=0.22, wspace=0.12)
                # figure(figsize=(3.5,4.5))
            
            
                close('all')
                rcParams['font.size'] = 8.
                rc('figure.subplot', left=0.2, right=0.93, bottom=0.12, top=0.95,
                                    hspace=0.25, wspace=0.12)


                Tlevs = arange(0,7,0.5)
                Ylabs = [0,500,1000,1500,2000]
                VTlevs = 1e-3*(arange(-20,20,2)+1)
                VTticks = 1e-3*(arange(-16,17,4))

                #figure(1, figsize=(3.5,6))
                figure(1, figsize=(3.5,6.))
                ascale = 1e-7
                aspace = 20

                Tlevs = arange(0,7,0.5)
                Ylabs = [0,500,1000,1500,2000]
                VTlevs = 1e-3*(arange(-20,20,2)+1)
                VTticks = 1e-3*(arange(-16,17,4))
                figure(1)
                ascale = 1e-7
                aspace = 20
            
                Tbar = T.mean(axis=2)
                subplot(3,1,2)
                cf=contourf(Y0, b.zc, Tflux_g_t / b.Lx, VTlevs, cmap=get_cmap('posneg'), extend='both')
                #clim([-0.02,0.02])
                #quiver(b.yc[::aspace] , b.zc[kr], vt[kr][:,::aspace], wt[kr][:,::aspace],
                #    angles='xy', scale_units='xy', scale=ascale)
                #c=contour(b.yc , b.zc,Tbar, Tlevs, colors='0.5')
                c=contour(Y00, b.zc, Tbar_T0, Tlevs, colors='0.5') 
                clabel(c,[0.5,], fmt='%1.1f')
                #title(r"$\bar{v}_g \bar{\theta}$ (Kms$^{-1}$) across %s" % mytit)
                title(r"$\bar{\mathbf{u}_g} \bar{\theta} \cdot \hat{\mathbf{n}}_{%s}$ (Kms$^{-1}$)" % sname)
                ylabel('Z (m)')
                xlim([0,2000e3])
                ylim([-2000,0])
                xticks(array(Ylabs)*1e3, [])
                Tbar = T.mean(axis=2)
                subplot(3,1,3)
                cf=contourf(Y0, b.zc, Tflux_e_t / b.Lx, VTlevs, cmap=get_cmap('posneg'), extend='both')
                #clim([-0.02,0.02])
                #quiver(b.yc[::aspace] , b.zc[kr], vt[kr][:,::aspace], wt[kr][:,::aspace],
                #    angles='xy', scale_units='xy', scale=ascale)
                #c=contour(b.yc , b.zc,Tbar, Tlevs, colors='0.5')
                c=contour(Y00, b.zc, Tbar_T0, Tlevs, colors='0.5') 
                clabel(c,[0.5,], fmt='%1.1f')
                #title(r"$\bar{v'\theta'}$ (Kms$^{-1}$) across %s" % mytit)
                title(r"$\bar{\mathbf{u}' \theta'} \cdot \hat{\mathbf{n}}_{%s}$ (Kms$^{-1}$)" % sname)
                ylabel('Z (m)')
                xlim([0,2000e3])
                ylim([-2000,0])
                xticks(array(Ylabs)*1e3, Ylabs)
                xlabel(r'$Y_{eq}$ (km)')
                #gca().fill_between(b.yc[r_[0,-1]], -array([2058.0,2058.0]),-array([b.H,b.H]),
                #    edgecolor='none', facecolor='0.7', alpha=0.3)       
                cb=colorbar(cf, cax=axes([0.05,0.04,0.9,0.01]), ticks=VTticks[::2], orientation='horizontal')
                #savefig('%s/heat_flux_%s-%s.pdf' % (figdir, avg, r))


                # vertical integrals

                tit = r'$\tau_0 =$ %4.3f N m$^{-2}$, %s contours, %s' % (b.tau0, avg, b.bathy)
                #figure(figsize=(3.5,2))
                subplot(3,1,1)
                #clf()
                plot(Y0, H_Ek /1e12, 'k',
                     Y0, H_SE /1e12, 'c',
                     Y0, H_TE /1e12, 'm',
                     Y0, (H_Ek + H_SE + H_TE) /1e12, '0.7'
                    )
                #grid()
                myrun = run
                if myrun=='bump':
                    myrun='ridge'
                if avg=='LAT':
                    xlabel('y (km)')
                    hname = r'\mathcal{H}'
                    mytit = 'Latitude Circles : %s' % myrun
                elif avg=='THETA':
                    xlabel(r'$\Theta$ ($^\circ$ C)')
                    hname = r'\mathcal{H}^\Theta'
                    mytit = r'$\Theta$ Contours : %s' % myrun
                elif avg=='PSI':
                    xlabel(r'$\Psi$ (Sv)')        
                    hname = r'\mathcal{H}^\Psi'
                    mytit = r'$\Psi$ Contours : %s' % myrun
                    if b.tau0==0.2:
                        xlim([-40,90])
                xlabel('')
                ylabel(r'$\mathcal{H}$ (TW)')
                xticks(array(Ylabs)*1e3, [])
                if b.tau0==0.2:
                    ylim([-120,120])
                legend([r'$%s_{Ek}$' % hname, r'$%s_{SE}$' % hname, r'$%s_{TE}$' % hname], 'upper left')
                #legend([r'$%s_{Ek}$' % hname, r'$%s_g$' % hname, r'$%s_{SE}$' % hname,
                #        r'$%s_{TE}$' % hname, r'$%s$' % hname], 'upper left')
                title('Heat Transport Across ' + mytit)
                #tight_layout()
                savefig('%s/heat_transport_%s-%s.pdf' % (figdir, avg, r))
        


#np.savez('./data/heat_transport.npz', HT=HT, tau=mytau)
#np.savez('./data/heat_transport_full_new.npz', HT=HT,
#        HT_coords=HT_coords, HT_eqlat=HT_eqlat, tau=mytau)

# streamwise coordinates
if False:
    X = linspace(0,b.Lx,b.Nx)/1e3
    Y = linspace(0,b.Ly,b.Ny)/1e3
    PSI = cumsum(b.integrate_vertical(U), axis=0) * b.dyc[0] / 1e6
    THETA = b.average_vertical(T)
    figure(figsize=(3.5,4.5))
    subplot(2,1,1)
    c=contour(X,Y, THETA, colors='k')
    clabel(c, fmt='%1.1f', fontsize=7)
    xlabel('x (km)'); ylabel('y (km)');
    title(r'$\Theta$ ($^\circ$ C)')
    subplot(2,1,2)
    c=contour(X,Y, PSI, colors='k')
    clabel(c, fmt='%1.1f', fontsize=7)
    xlabel('x (km)'); ylabel('y (km)');
    title(r'$\Psi$ (Sv)')
    tight_layout()
    savefig('figures/paper/streamwise_coords.pdf')



