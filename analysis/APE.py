from pylab import *
import bump_channel
#import mycolors
import os

base_dir = os.path.join(os.environ['D'],'DATASTORE.RPA','projects','drag_strat')

b_bump = bump_channel.BumpChannel(base_dir, 'taux2000_rb0110_bump')
b_flat = bump_channel.BumpChannel(base_dir, 'taux2000_rb0110_flat')

for B in [b_bump]: #,b_flat]:
    c = B.get_APEcomputer()
    Ep = c.calc_PE()
    Eb = c.calc_BPE()
    Ea = c.calc_APE()
    
    # Stu's QG approximation
    bmean = c.b.mean(axis=1).mean(axis=1)
    N2 = B.ddz_cgrid_centered(bmean)
    Ea_QG = sum(c.vol * 0.5 * (c.b-bmean[:,newaxis,newaxis])**2 / N2[:,newaxis,newaxis])
    
    bbar = tile(c.b.mean(axis=2)[:,:,newaxis],[1,1,B.Nx])
    bbar_z = B.ddz_cgrid_centered(bbar)
    bwave = c.b - bbar
    
    print('%040s: %4.2f PJ' % ('Total Potential Energy',Ep/1e15))
    print('%040s: %4.2f PJ' % ('Background Potential Energy',Eb/1e15))
    print('%040s: %4.2f PJ' % ('Available Potential Energy',Ea/1e15))
    print('%040s: %4.2f PJ' % ('Available Potential Energy (QG)',Ea_QG/1e15))

    BstarC = zeros((B.Nz,B.Ny))
    for j in arange(B.Ny):
        zf_vol = hstack([0,cumsum(c.vol[:,j].sum(axis=1))])
        if zf_vol[-1] > 0:
            BstarF = c.calc_Bstar_f(c.b[:,j],c.vol[:,j],zf_vol)
            BstarC[:,j] = 0.5 * (BstarF[1:] + BstarF[:-1])
    
    Ea_wave = sum(c.vol * B.zc[:,newaxis,newaxis] * -(c.b-BstarC[:,:,newaxis]))
    Ea_wave_QG = sum(c.vol * 0.5 * (bwave**2) / bbar_z)

    print('%040s: %4.2f PJ' % ('APE in standing wave',Ea_wave/1e15))
    print('%040s: %4.2f PJ' % ('APE in standing wave (QG)',Ea_wave_QG/1e15))
    