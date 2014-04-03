from pylab import *
import os
import bump_channel
import MITgcmdata.poisson
import mycolors

rcParams['figure.figsize'] = [6, 3.5]
rcParams['legend.fontsize'] = 7
rcParams['font.size'] = 8

base_dir = os.path.join(os.environ['D'],'DATASTORE.RPA','projects','drag_strat')
figdir = './figures/runs_HT'

#mytau = array([0.0125, 0.0250, 0.05, 0.1, 0.2, 0.4, 0.8])
myrb = array([ 0.00014,  0.00028,  0.00055,  0.0011 ,  0.0022 ,  0.0044 ,  0.0088 ])
tau = 0.1
runs = ['flat', 'bump']
pathtypes = ['LAT', 'PSI', 'THETA']
pathnames = ['Lat', r'$\Psi$', r'$\Theta$']
components = ['Mean','Standing','Transient']

# spacing in streamline space
Nt = 100

# heat transport arrays
HT = zeros((len(runs), len(pathtypes), len(components), len(myrb), Nt))
HT_coords = zeros((len(runs), len(pathtypes), len(myrb), Nt))
HT_eqlat = zeros((len(runs), len(pathtypes), len(myrb), Nt))

#for nrun in arange(len(runs)):
#for nrun in [0,]:
for nrun in [1,]:
    run = runs[nrun]
    for nrb in arange(len(myrb)):
#    for nrb in arange(len(mytau)):
#    for nrb in [4,]:

        close('all')

        rb = myrb[nrb]
        r = 'taux%04d_rb%04d_%s' % (tau*10000, ceil(rb*1e5), run)
        b = bump_channel.BumpChannel(base_dir, r)
        sol = MITgcmdata.poisson.solver(b) # conjugate gradient solver

        T = b.rdmds('Ttave')
        Tw = b.cgrid_to_ugrid(T) # on the U grid
        Ts = b.cgrid_to_vgrid(T) # on the V grid
        V = b.rdmds('Vveltave')
        U = b.rdmds('Uveltave')
        P = b.rdmds('Phhytave')

        EKE = 0.5*b.average_vertical((b.rdmds('UUtave')-U**2+b.rdmds('VVtave')-V**2))

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
        UTtave = b.integrate_vertical(U*Tw)
        VTtave = b.integrate_vertical(V*Ts)
        divTflux = b.ddx_ugrid_to_cgrid(UTtave) + b.ddy_vgrid_to_cgrid(VTtave)
        #UTtave_T0 = b.integrate_vertical(U*(Tw-THETAw))
        #VTtave_T0 = b.integrate_vertical(V*(Ts-THETAs))
        #divTflux_T0 = b.ddx_ugrid_to_cgrid(UTtave_T0) + b.ddy_vgrid_to_cgrid(VTtave_T0)
        # ageostrophic flux
        UTa = b.integrate_vertical(Ua*Tw)
        VTa = b.integrate_vertical(Va*Ts)
        divTflux_a = b.ddx_ugrid_to_cgrid(UTa) + b.ddy_vgrid_to_cgrid(VTa)
        #UTa_T0 = b.integrate_vertical(Ua*(Tw-THETAw))
        #VTa_T0 = b.integrate_vertical(Va*(Ts-THETAs))
        #divTflux_a_T0 = b.ddx_ugrid_to_cgrid(UTa_T0) + b.ddy_vgrid_to_cgrid(VTa_T0)
        # geostrophic flux
        UTg = b.integrate_vertical(Ug*Tw)
        VTg = b.integrate_vertical(Vg*Ts)
        divTflux_g = b.ddx_ugrid_to_cgrid(UTg) + b.ddy_vgrid_to_cgrid(VTg)
        #UTg_T0 = b.integrate_vertical(Ug*(Tw-THETAw))
        #VTg_T0 = b.integrate_vertical(Vg*(Ts-THETAs))
        #divTflux_g_T0 = b.ddx_ugrid_to_cgrid(UTg_T0) + b.ddy_vgrid_to_cgrid(VTg_T0)
        #divTflux_g_BT = ( b.ddx_ugrid_to_cgrid(b.integrate_vertical(Ug) * THETAw) +
        #                  b.ddy_vgrid_to_cgrid(b.integrate_vertical(Vg) * THETAs) )
        # eddy flux
        UTe = b.integrate_vertical(UT - U*Tw)
        VTe = b.integrate_vertical(VT - V*Ts)
        divTflux_e = b.ddx_ugrid_to_cgrid(UTe) + b.ddy_vgrid_to_cgrid(VTe)

        # surface eddy flux
        UTe0, VTe0 = (UT - U*Tw)[0], (VT - V*Ts)[0]
        divTflux_e0 = b.ddx_ugrid_to_cgrid(UTe0) + b.ddy_vgrid_to_cgrid(VTe0)
        
        # for K
        kbot = 30
        Tbot = 0
        D = -2 * b.integrate_vertical( b.zc[:,newaxis,newaxis] * (T - Tbot), krange=r_[0:kbot]) /  (
            b.integrate_vertical( T.mean(axis=2).mean(axis=1)[:,newaxis]-Tbot, krange=r_[0:kbot]) )
        Ttc = zeros((b.Ny,b.Nx))
        for j in arange(b.Ny):
            for i in arange(b.Nx):
                Ttc[j,i] = mean(interp(linspace(0,D[j,i],100), -b.zc, T[:,j,i]))
        Ttc = ma.masked_array(Ttc, T[0].mask)
                
        #Zgrid = tile(b.zc[:,newaxis,newaxis],(1,b.Ny,b.Nx))
        #Ttc = b.integrate_vertical(ma.masked_array(T, Zgrid < -D[newaxis,:,:])) / D

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
                    pass
                    #extra = tmean.max() * log2(tau/mytau[0])/7.5
                t = linspace(tmean.min()-extra, tmean.max()+extra, Nt)
                #t = linspace(tmean.min()-extra, tmean.max(), Nt)
            
   
            # output arrays
            Tflux_a_t = zeros(Nt)
            Tflux_g_t = zeros(Nt)
            Tflux_e_t = zeros(Nt)

            # divergence method
            tbnds = hstack([-inf, 0.5*(t[1:]+t[:-1]), inf])
            # special value at which to do the integral as a function of T
            n0 = 40 # near 0.9 deg. isotherm
            for n in arange(Nt):
                idx = T0<=tbnds[n+1]
                area = b.rac[idx].sum()
                Tflux_a_t[n] = (b.rac*divTflux_a)[idx].sum() 
                Tflux_g_t[n] = (b.rac*divTflux_g)[idx].sum()
                Tflux_e_t[n] = (b.rac*divTflux_e)[idx].sum()
                if n==n0:
                    Tflux_a_T0_x = zeros(b.Nx)
                    Tflux_g_T0_x = zeros(b.Nx)
                    Tflux_e_T0_x = zeros(b.Nx)                    
                    for i in arange(b.Nx):
                        midx = idx & (b.xc[newaxis,:]<=b.xc[i])
                        Tflux_a_T0_x[i] = (b.rac*divTflux_a)[midx].sum() 
                        Tflux_g_T0_x[i] = (b.rac*divTflux_g)[midx].sum()
                        Tflux_e_T0_x[i] = (b.rac*divTflux_e)[midx].sum()                     


            # normal flux method
            if True:
                # for some reason, this doesn't always work
                sol = MITgcmdata.poisson.solver(b)
                sol.solve(divTflux.copy())
                UT_div, VT_div = sol.get_vectors()
                sol = MITgcmdata.poisson.solver(b)
                sol.solve(divTflux_g.copy())
                UTg_div, VTg_div = sol.get_vectors()
                sol = MITgcmdata.poisson.solver(b)
                sol.solve(divTflux_e.copy())
                UTe_div, VTe_div = sol.get_vectors()
                #sol = MITgcmdata.poisson.solver(b)
                #sol.solve(divTflux_a.copy())
                #UTaa_div, VTaa_div = sol.get_vectors()
                #sol = MITgcmdata.poisson.solver(b)
                #sol.solve(divTflux_e0.copy())
                #UTe0_div, VTe0_div = sol.get_vectors()



                UTa_div, VTa_div = UT_div-UTg_div, VT_div-VTg_div


                # normal coordinate
                T0x = b.ddx_cgrid_centered(T0)
                T0y = b.ddy_cgrid_centered(T0)
                gradT = (T0x**2 + T0y**2)**(0.5)
                # project everything normal
                Tflux_a = b.ugrid_to_cgrid(UTa) * T0x / gradT + b.vgrid_to_cgrid(VTa) * T0y / gradT
                Tflux_g = b.ugrid_to_cgrid(UTg) * T0x / gradT + b.vgrid_to_cgrid(VTg) * T0y / gradT
                Tflux_e = b.ugrid_to_cgrid(UTe) * T0x / gradT + b.vgrid_to_cgrid(VTe) * T0y / gradT
                # divergent component only
                Tflux_a_div = b.ugrid_to_cgrid(UTa_div) * T0x / gradT + b.vgrid_to_cgrid(VTa_div) * T0y / gradT
                Tflux_g_div = b.ugrid_to_cgrid(UTg_div) * T0x / gradT + b.vgrid_to_cgrid(VTg_div) * T0y / gradT
                Tflux_e_div = b.ugrid_to_cgrid(UTe_div) * T0x / gradT + b.vgrid_to_cgrid(VTe_div) * T0y / gradT
                # surface component
                #T00x, T00y = b.ddx_cgrid_centered(T[0]), b.ddy_cgrid_centered(T[0])
                #gradT00 = (T00x**2 + T00y**2)**(0.5)
                #Tflux_e0_div = b.ugrid_to_cgrid(UTe0_div) * T00x / gradT00 + b.vgrid_to_cgrid(VTe0_div) * T00y / gradT00
        
                # gradient just over thermocline
                gradTtc = (b.ddx_cgrid_centered(Ttc)**2 + b.ddy_cgrid_centered(Ttc)**2)**0.5

                # diffusivities
                Kmask = log10(gradT**-1) > 7.5 # maks weak gradient areas
                K = - ma.masked_array(Tflux_e / gradT / D, Kmask)
                Kdiv = - ma.masked_array(Tflux_e_div / gradT / D, Kmask)
                Kdiv_tc = - ma.masked_array(Tflux_e_div / gradTtc / D, Kmask)
                Kdiv_H = - ma.masked_array(Tflux_e_div / gradT / b.depth, Kmask)
                #K0 = - ma.masked_array(Tflux_e0_div/gradT00)
                

                rc('figure.subplot', left=0.04, right=0.92, bottom=0.1, top=0.78, wspace=0.2, hspace=0.4)
                figure(1,figsize=(6.5,2.2))
                clf()       
                mylims = array([-1,1]) * 2 * b.tau0  * 8 / (-b.f0 * b.rho0)
                fluxlevs = (arange(-30,30,5)+2.5) # * rb/1.1e-3
                fluxticks = arange(-25,26,5) # * rb/1.1e-3
                adens = 14 # arrow density
                ascale = 600 # * rb/1.1e-3
                
                subplot(131)
                #pcolormesh(Tflux_a_div, cmap=get_cmap('posneg'), rasterized=True); clim(mylims)
                cf=contourf(b.xc, b.yc, Tflux_a_div, fluxlevs, cmap=get_cmap('posneg'), extend='both', rasterized=True);
                c=contour(b.xc, b.yc, T0, colors='0.5')
                quiver(b.xc[(adens/2)::adens],b.yc[(adens/2)::adens],
                    UTa_div[(adens/2)::adens,(adens/2)::adens],VTa_div[(adens/2)::adens,(adens/2)::adens], scale=ascale)               
                clabel(c, fmt='%1.1f', fontsize=7)
                #title(r'$\int dz (\bar {\bf v} \bar \theta)_a^{div} \cdot \hat n$')
                title(r'$(\bar {\bf v} \bar \theta)_a^{div}$')
                xticks([]); xlabel('x'); yticks([]); ylabel('y')
                
                subplot(132)
                #pcolormesh(Tflux_g_div, cmap=get_cmap('posneg'), rasterized=True); clim(mylims)
                cf=contourf(b.xc, b.yc, Tflux_g_div, fluxlevs, cmap=get_cmap('posneg'), extend='both', rasterized=True);
                c=contour(b.xc, b.yc, T0, colors='0.5')
                quiver(b.xc[(adens/2)::adens],b.yc[(adens/2)::adens],
                    UTg_div[(adens/2)::adens,(adens/2)::adens],VTg_div[(adens/2)::adens,(adens/2)::adens], scale=ascale)
                clabel(c, fmt='%1.1f', fontsize=7)
                #title(r'$\int dz (\bar {\bf v} \bar \theta)_g^{div} \cdot \hat n$')
                title(r'$(\bar {\bf v} \bar \theta)_g^{div}$')
                xticks([]); xlabel('x'); yticks([]); ylabel('y')

                subplot(133)
                #pcolormesh(Tflux_e_div, cmap=get_cmap('posneg'), rasterized=True); clim(mylims)
                cf=contourf(b.xc, b.yc, Tflux_e_div, fluxlevs, cmap=get_cmap('posneg'), extend='both', rasterized=True);
                c=contour(b.xc, b.yc, T0, colors='0.5')
                qu=quiver(b.xc[(adens/2)::adens],b.yc[(adens/2)::adens],
                    UTe_div[(adens/2)::adens,(adens/2)::adens],VTe_div[(adens/2)::adens,(adens/2)::adens], scale=ascale)
                clabel(c, fmt='%1.1f', fontsize=7)
                #title(r"$\int dz (\bar{ {\bf v'} \theta'})^{div} \cdot \hat n$")
                title(r"$(\bar{ {\bf v'} \theta'})^{div}$")
                xticks([]); xlabel('x'); yticks([]); ylabel('y')

                cb=colorbar(cf, ticks=fluxticks, cax=axes([0.93,0.1,0.01,0.68]))
                cl = getp(cb.ax, 'ymajorticklabels') 
                setp(cl, fontsize=7)
                qk = quiverkey(qu, 0.06, 0.15, 30, r'30 K m$^2$ s$^{-1}$',
                               labelpos='E', coordinates='figure', fontproperties={'size':7})
                
                gcf().text(0.5,0.92,r'Vertically Integrated Divergent Heat Fluxes (K m$^2$ s$^{-1}$)', ha='center')
                
                #tight_layout()
                savefig('%s/normal_flux_div_%s-%s.pdf' % (figdir, avg, r))

                # just two figures
                fluxlevs = (arange(-25,25,5)+2.5) # * rb/1.1e-3
                fluxticks = arange(-20,21,5) # * rb/1.1e-3
                mytick = [500,1000,1500] 
                rc('figure.subplot', left=0.09, right=0.88, bottom=0.15, top=0.82, wspace=0.12, hspace=0.12)
                figure(33,figsize=(7,3.0))
                clf()       
                
                subplot(121)
                #axes([0.09,0.15,0.25,0.65])
                cf=contourf(b.xc, b.yc, Tflux_a_div + Tflux_g_div, fluxlevs,
                    cmap=get_cmap('posneg'), extend='both', rasterized=True);
                c=contour(b.xc, b.yc, T0, colors='0.5')
                quiver(b.xc[(adens/2)::adens],b.yc[(adens/2)::adens],
                    (UTa_div+UTg_div)[(adens/2)::adens,(adens/2)::adens],
                    (VTa_div+VTg_div)[(adens/2)::adens,(adens/2)::adens], scale=ascale)               
                clabel(c, fmt='%1.1f', fontsize=7)
                #title(r'$\int dz (\bar {\bf v} \bar \theta)_a^{div} \cdot \hat n$')
                title(r'$(\bar {\bf v} \bar \theta)^{div}$')
                #xticks([]); xlabel('x'); yticks([]); ylabel('y')
                xticks(array(mytick)*1000, mytick); yticks(array(mytick)*1000, mytick)
                xlabel('x (km)'); ylabel('y (km)')
                
                subplot(122)
                #axes([0.35,0.15,0.25,0.65])
                cf=contourf(b.xc, b.yc, Tflux_e_div, fluxlevs, cmap=get_cmap('posneg'), extend='both', rasterized=True);
                c=contour(b.xc, b.yc, T0, colors='0.5')
                qu=quiver(b.xc[(adens/2)::adens],b.yc[(adens/2)::adens],
                    UTe_div[(adens/2)::adens,(adens/2)::adens],VTe_div[(adens/2)::adens,(adens/2)::adens], scale=ascale)
                clabel(c, fmt='%1.1f', fontsize=7)
                #title(r"$\int dz (\bar{ {\bf v'} \theta'})^{div} \cdot \hat n$")
                title(r"$(\bar{ {\bf v'} \theta'})^{div}$")
                #xticks([]); xlabel('x'); yticks([]); ylabel('y')
                xticks(array(mytick)*1000, mytick); yticks(array(mytick)*1000, [])
                xlabel('x (km)');
                
                #subplot(133)
                # Qfac = b.rho0 * 3994  / 1e12
                # axes([0.67,0.15,0.25,0.65])
                # plot(b.xc, Qfac*(Tflux_a_T0_x + Tflux_g_T0_x), 'k')
                # plot(b.xc, Qfac*(Tflux_e_T0_x), 'c')
                # plot(b.xc, Qfac*(Tflux_a_T0_x + Tflux_g_T0_x + Tflux_e_T0_x), '0.7')
                # xticks(array(mytick)*1000, mytick);
                # title('Cumulative Integral in x')
                # legend(['mean','eddy'], loc='upper left')
                # grid()
                # xlabel('x (km)');
                # ylabel('Heat Transport (TW)')
                # gca().yaxis.tick_right()
                # gca().yaxis.set_label_position('right')

                cb=colorbar(cf, ticks=fluxticks, cax=axes([0.91,0.15,0.01,0.67]))
                cl = getp(cb.ax, 'ymajorticklabels') 
                setp(cl, fontsize=7)
                qk = quiverkey(qu, 0.12, 0.2, 30, r'30 K m$^2$ s$^{-1}$',
                               labelpos='E', coordinates='figure', fontproperties={'size':7})
                
                gcf().text(0.5,0.91,r'Vertically Integrated Divergent Heat Fluxes (K m$^2$ s$^{-1}$)', ha='center')
                
                #tight_layout()
                savefig('%s/normal_flux_2_div_%s-%s.pdf' % (figdir, avg, r))
                
                #Klevs = arange(-4000,4000,500)+250
                #Kticks = arange(-3500,3501,500)
                Klevs = arange(0,5000,500) # * rb/1.1e-3
                Kticks = arange(0,5000,1000) # * rb/1.1e-3
                EKElevs = arange(0,31,2)
                qscale = ascale
                
                #figure(2,figsize=(6.5,2.2))    
                figure(2, figsize=(3.5,2.2))
                #Klevs = arange(-7000,7001,2000)

                #pcolormesh(b.xc, b.yc, Kdiv, cmap=get_cmap('posneg'), rasterized=True); clim(array([-1,1])*8e3)
                pc=contourf(b.xc, b.yc, Kdiv, Klevs, cmap=get_cmap('jet'), extend='both', rasterized=True);
                c=contour(b.xc, b.yc, T0, colors='0.5')
                clabel(c, fmt='%1.1f', fontsize=7)
                quiver(b.xc[(adens/2)::adens],b.yc[(adens/2)::adens],
                    UTe_div[(adens/2)::adens,(adens/2)::adens],VTe_div[(adens/2)::adens,(adens/2)::adens], scale=qscale)
                title(r"$K_\perp^{div}$ (m$^2$ s$^{-1}$) Thermoclice Average")
                xticks([]); xlabel('x'); yticks([]); ylabel('y')
                tight_layout()
                cb=colorbar(pc, ticks=Kticks)
                #cb=colorbar(pc, cax=axes([0.91,0.1,0.01,0.8]))
                cl = getp(cb.ax, 'ymajorticklabels') 
                setp(cl, fontsize=7)
                
                savefig('%s/Kdiv_tc_%s-%s.pdf' % (figdir, avg, r))

                figure(3, figsize=(3.5,2.2))
                #Klevs = arange(-7000,7001,2000)

                #pcolormesh(b.xc, b.yc, Kdiv, cmap=get_cmap('posneg'), rasterized=True); clim(array([-1,1])*8e3)
                pc=contourf(b.xc, b.yc, Kdiv_H, Klevs, cmap=get_cmap('jet'), extend='both', rasterized=True);
                c=contour(b.xc, b.yc, T0, colors='0.5')
                clabel(c, fmt='%1.1f', fontsize=7)
                quiver(b.xc[(adens/2)::adens],b.yc[(adens/2)::adens],
                    UTe_div[(adens/2)::adens,(adens/2)::adens],VTe_div[(adens/2)::adens,(adens/2)::adens], scale=qscale)
                title(r"$K_\perp^{div}$ (m$^2$ s$^{-1}$) Depth Average")
                xticks([]); xlabel('x'); yticks([]); ylabel('y')
                #cb=colorbar(pc, cax=axes([0.91,0.1,0.01,0.8]))
                tight_layout()
                cb=colorbar(pc, ticks=Kticks)
                cl = getp(cb.ax, 'ymajorticklabels') 
                setp(cl, fontsize=7)
                
                savefig('%s/Kdiv_depth_%s-%s.pdf' % (figdir, avg, r))                

                figure(4, figsize=(3.5,2.2))
                # #pcolormesh(b.xc, b.yc, Kdiv, cmap=get_cmap('posneg'), rasterized=True); clim(array([-1,1])*8e3)
                pc=contourf(b.xc, b.yc, EKE*1000, EKElevs, cmap=get_cmap('jet'), extend='both', rasterized=True);
                c=contour(b.xc, b.yc, T[0], colors='0.5')
                clabel(c, fmt='%1.1f', fontsize=7)
                quiver(b.xc[(adens/2)::adens],b.yc[(adens/2)::adens],
                    UTe_div[(adens/2)::adens,(adens/2)::adens],VTe_div[(adens/2)::adens,(adens/2)::adens], scale=500)
                title(r"EKE (cm$^2$ s$^{-2}$) (Vertical Avg)")
                xticks([]); xlabel('x'); yticks([]); ylabel('y')
                # #cb=colorbar(pc, cax=axes([0.91,0.1,0.01,0.8]))
                tight_layout()
                cb=colorbar(pc) #, ticks=(2*Kticks))
                cl = getp(cb.ax, 'ymajorticklabels') 
                setp(cl, fontsize=7)
                #             
                savefig('%s/EKE_%s-%s.pdf' % (figdir, avg, r))                
        
            tit = r'$\tau_0 =$ %4.3f N m$^{-2}$, %s contours, %s' % (b.tau0, avg, b.bathy)
            figure(figsize=(3.5,2))
            clf()
            Qfac = b.rho0 * 3994  / 1e12
            plot(t, Qfac*Tflux_a_t, 'k',
                 t, Qfac*(Tflux_g_t + Tflux_e_t), 'c',
                 t, Qfac*Tflux_g_t, 'm',
                 t, Qfac*Tflux_e_t, 'y',
                 t, Qfac*(Tflux_g_t+Tflux_a_t+Tflux_e_t), '0.7'
                )
            grid()
            myrun = run
            if myrun=='bump':
                myrun='ridge'
            if avg=='LAT':
                xlabel('y (km)')
                hname = r'\mathcal{H}'
                mytit = 'Latitude Circles - %s' % myrun
            elif avg=='THETA':
                xlabel(r'$\Theta$ ($^\circ$ C)')
                hname = r'\mathcal{H}^\Theta'
                mytit = r'$\Theta$ Contours - %s' % myrun
            elif avg=='PSI':
                xlabel(r'$\Psi$ (Sv)')        
                hname = r'\mathcal{H}^\Psi'
                mytit = r'$\Psi$ Contours - %s' % myrun
                if b.tau0==0.2:
                    xlim([-40,90])
            ylabel('Heat Transport (TW)')
            if b.tau0==0.2:
                ylim([-120,120])
            legend([r'$%s_{Ek}$' % hname, r'$%s_g$' % hname, r'$%s_{SE}$' % hname,
                    r'$%s_{TE}$' % hname, r'$%s$' % hname], 'upper left')
            title(mytit)
            tight_layout()
            savefig('%s/heat_transport_%s-%s.pdf' % (figdir, avg, r))
        
            # calculate max ht
            j = argmax(Tflux_a_t)
            print 'Max j: %g' % j
            # HT = zeros((len(runs), len(pathtypes), len(components), len(mytau)))
            HT[nrun, npath, 0, nrb] = Qfac*Tflux_a_t
            HT[nrun, npath, 1, nrb] = Qfac*Tflux_g_t
            HT[nrun, npath, 2, nrb] = Qfac*Tflux_e_t
            HT_coords[nrun, npath, nrb] = t
            
            # calc equivalent y
            HT_eqlat[nrun, npath, nrb] = interp(t,T0.mean(axis=1),b.yc)

#np.savez('./data/heat_transport.npz', HT=HT, tau=mytau)
#np.savez('./data/heat_transport_full.npz', HT=HT,
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



