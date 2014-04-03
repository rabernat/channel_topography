from pylab import *
import bump_channel
#import mycolors
import os

base_dir = os.path.join(os.environ['D'],'DATASTORE.RPA','projects','drag_strat')

tau = [0, 0.0125, 0.0250, 0.05, 0.1, 0.2, 0.4, 0.8]
N = len(tau)

res = dict()
for k in ['D', 'D05', 'APE','PEtoKE','EPEtoEKE','PEtoIE','KEdiss','EKEdiss','EKEbot','EKE']:
    res[k] = zeros(N)

for n in arange(N):
    tau0 = tau[n]
    if tau0==0:
        suf = 'flatshort'
    else:
        suf = 'flat'
        
    b = bump_channel.BumpChannel(base_dir, 'taux%04d_rb0110_%s' % (tau0*10000, suf))
    b.load_default_fields()
    U2 = b.rdmds('UUtave')
    V2 = b.rdmds('VVtave')
    EKE = 0.5 * (U2 - b.U**2 + V2 - b.V**2).mean(axis=2)
    
    # thermocline depth
    D = -2 * b.integrate_vertical( b.zc[:,newaxis]
        * (b.Tbar - b.Tbar[-1])) /  b.integrate_vertical(b.Tbar)
    res['D'][n] = D[-2]
    
    # depth of 0.5 degree isotherm
    res['D05'][n] = -interp(0.5, b.Tbar[::-1,-2], b.zc[::-1])

    # available potential energy
    res['APE'][n] = b.get_APEcomputer().calc_APE()
    
    # energy conversion from potential to kinetic (W m^2)
    PEfac = b.rho0 * b.g * b.tAlpha 
    res['PEtoKE'][n] = PEfac * mean(b.average_vertical(b.WTbar))
    res['EPEtoEKE'][n] = PEfac * mean(b.average_vertical(b.WpTpbar))
    res['PEtoIE'][n] = PEfac * mean(b.average_vertical(b.Tdifbar))
    res['KEdiss'][n] = b.rho0 * b.rb * mean(U2[-1] + V2[-1]) / b.H
    res['EKEdiss'][n] = b.rho0 * b.rb * 2 * mean(EKE[-1]) / b.H
    
    # mean EKE
    res['EKEbot'][n] = mean(EKE[-1])
    res['EKE'][n] = mean(b.average_vertical(EKE))
    
    # contour(b.yc /1e6, b.zc, Tbar, arange(0,8,0.5))
    # contour(b.yc /1e6, b.zc, b_wind.Tbar, arange(0,8,0.5), linestyles=':')
    # xlabel('Y (km)')
    # ylabel('Z (m)')
    # 
    # figure()
    # pcolormesh(b.yc/1e6, b.zc, ma.masked_invalid(log10(c.mean(axis=2))), cmap=get_cmap('CMRmap_r'))
    # colorbar()
    # xlabel('Y (km)')
    # ylabel('Z (m)')
    # title(r'log$_{10}$(convective frequency)')

rcParams['font.size'] = 8
rcParams['legend.fontsize'] = 8
figure(figsize=(3,4))
subplot(211)
plot(tau,res['APE'], 'k.-')
title('APE (J)')
xlabel(r'$\tau_0$ (N m$^{-2}$)')
xlim([-0.05,0.85])
grid()
subplot(212)
plot(tau,res['D'], 'k.-')
plot(tau,res['D05'], 'k*--')
legend(['moment','0.5 isotherm'], loc='lower right')
title('N.B. Thermocline Depth (m)')
xlabel(r'$\tau_0$ (N m$^{-2}$)')
xlim([-0.05,0.85])
grid()
tight_layout()
savefig('figures/nowind/APE_h_vs_tau0.pdf')






