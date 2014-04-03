from pylab import *
import bump_channel
import os
import mycolors
import MITgcmdata.poisson

# plot settings
figdir = './figures/'

# the list of runs
mytau = array([0.0125, 0.0250, 0.05, 0.1, 0.2, 0.4, 0.8])
runs = ['flat', 'bump']

base_dir = os.path.join(os.environ['D'],'DATASTORE.RPA','projects','drag_strat')

Nx = 400
N = len(mytau)

# length of streamwise coord
#Ns = 800
# use same streamwise coordinate for all runs
# some will have longer streamlines than others
Ns = 800
S0 = linspace(0,3.7e6,Ns)

def avg_cross_stream(q,s,S0):
    # s should be a masked array
    Ns = len(S0)
    qsum = zeros(Ns)
    area = zeros(Ns)
    for n in arange(Ns):
        fullmask = s.mask | (s>S0[n])
        area[n] = (1-fullmask).sum()
        qsum[n] = sum(ma.masked_array(q,fullmask))
    #subplot(121); plot(area); subplot(122); plot(qsum)
    Dqsum = diff(qsum)
    Darea = diff(area)
    qS = Dqsum/Darea
    # interpolate nans
    mask = isnan(qS)
    qS[mask] = interp(flatnonzero(mask), flatnonzero(~mask), qS[~mask])
    return ma.masked_array(hstack([qS, 0]), S0>s.max())


figure(figsize=(6.5,2.8))

n=-1
#for nrun in [1,]:
for nrun in [0,1]:
    run = runs[nrun]
    #for ntau in arange(N):
    for ntau in [4]:
        n+=1

        tau = mytau[ntau]
        r = 'taux%04d_rb0110_%s' % (tau*10000, run)
        try:
            b.reinit(r)
        except NameError:
            b = bump_channel.BumpChannel(base_dir, r)
        sol = MITgcmdata.poisson.solver(b)
        
        T = b.rdmds('Ttave')
        U = b.rdmds('Uveltave')
        V = b.rdmds('Vveltave')
        UT = b.rdmds('UTtave')
        VT = b.rdmds('VTtave')
        Vp2 =   b.rdmds('UUtave')-U**2 + b.rdmds('VVtave') - V**2
        Tp2 = b.rdmds('TTtave')-T**2
        X,Y = meshgrid(b.xc,b.yc)
        
        UpTp = UT - U*b.cgrid_to_ugrid(T)
        VpTp = VT - V*b.cgrid_to_vgrid(T)
        
        div_eddy_flux = b.integrate_vertical(
                     b.ddx_ugrid_to_cgrid(UpTp) +
                     b.ddy_vgrid_to_cgrid(VpTp) )

        #while (err2**0.5 > 2):
        sol = MITgcmdata.poisson.solver(b)
        sol.solve(div_eddy_flux.copy())
        UpTp_div, VpTp_div = sol.get_vectors()
    
        err2 = mean(( (b.ddx_ugrid_to_cgrid(UpTp_div)
               +b.ddy_vgrid_to_cgrid(VpTp_div))[2:]
               -div_eddy_flux[2:] )**2) / mean(div_eddy_flux**2)
           
        print 'Normalized error: %g' % err2**0.5
        print 'Solvers error %g' % sol.mynorm
        print 'Congergence ratio %g' % (sol.residuals[-1]/sol.residuals[0])**(1.0/len(sol.residuals))

        T0 = b.average_vertical(T)
        T0x = b.ddx_cgrid_centered(T0)
        T0y = b.ddy_cgrid_centered(T0)
        gradT = (T0x**2 + T0y**2)**(0.5)
        
        # project everything normal
        # divergent component only
        Tflux_div = b.ugrid_to_cgrid(UpTp_div) * T0x / gradT + b.vgrid_to_cgrid(VpTp_div) * T0y / gradT
        
        # mixing efficiency parameter
        alph = -Tflux_div / b.integrate_vertical(Vp2 * Tp2)
        
        subplot(1,2,n+1)
        pcolormesh(b.xc/1e6, b.yc/1e3, log10(ma.masked_less_equal( alph, 0)), rasterized=True)
        clim([-2,1])
        colorbar()
        contour(b.xc/1e6, b.yc/1e3, T0, 8, colors='k')
        xlabel('x (km)'); ylabel('y (km)')
        title(runs[nrun])

t = gcf().text(0.5,0.92,r"log$_{10}(-\overline{v'\theta'} / (v'_{rms} \theta'_{rms} ) )$", ha='center')
        
