from pylab import *
import bump_channel
import os
import mycolors
import MITgcmdata.poisson
import MITgcmdata.poisson2

# plot settings
figdir = './figures/'
base_dir = os.path.join(os.environ['D'],'DATASTORE.RPA','projects','drag_strat')

tau=0.2
r = 'taux2000_rb0110_bump'
b = bump_channel.BumpChannel(base_dir, r)
sol = MITgcmdata.poisson.solver(b)
sol2 = MITgcmdata.poisson2.solver(b)

T = b.rdmds('Ttave')
U = b.rdmds('Uveltave')
V = b.rdmds('Vveltave')
UT = b.rdmds('UTtave')
VT = b.rdmds('VTtave')
X,Y = meshgrid(b.xc,b.yc)

UpTp = UT - U*b.cgrid_to_ugrid(T)
VpTp = VT - V*b.cgrid_to_vgrid(T)

# how to calculate divergence according to MITgcm source code
#          uTrans(i,j,bi,bj) = dyG(i,j,bi,bj)*drF(k)
# &                          *uvel(i,j,ks,bi,bj)*maskInW(i,j,bi,bj)
#          vTrans(i,j,bi,bj) = dxG(i,j,bi,bj)*drF(k)
# &                          *vvel(i,j,ks,bi,bj)*maskInS(i,j,bi,bj)
# C-    calculate RHS = rAc*Div(uVel,vVel):
#               b2d(i,j,bi,bj) = (
#      &                    ( uTrans(i+1,j,bi,bj) - uTrans(i,j,bi,bj) )
#      &                  + ( vTrans(i,j+1,bi,bj) - vTrans(i,j,bi,bj) )
#      &                         )*maskInC(i,j,bi,bj)
     
fluxX = sum( UpTp * b.dyg[newaxis,newaxis,:] * b.dzf[:,newaxis,newaxis]  * b.hFacW, axis=0 )
fluxY = sum( VpTp * b.dxg[newaxis,newaxis,:] * b.dzf[:,newaxis,newaxis]  * b.hFacS, axis=0 )

div = -fluxX + roll(fluxX,-1,axis=1) - fluxY + roll(fluxY,-1,axis=0)

sol2.solve(div)
u,v = sol2.get_vectors()

# div_eddy_flux = b.integrate_vertical(
#              b.ddx_ugrid_to_cgrid(UpTp) +
#              b.ddy_vgrid_to_cgrid(VpTp) )
# 
# #while (err2**0.5 > 2):
# sol = MITgcmdata.poisson.solver(b)
# sol.solve(div_eddy_flux.copy())
# UpTp_div, VpTp_div = sol.get_vectors()
# 
# err2 = mean(( (b.ddx_ugrid_to_cgrid(UpTp_div)
#        +b.ddy_vgrid_to_cgrid(VpTp_div))[2:]
#        -div_eddy_flux[2:] )**2) / mean(div_eddy_flux**2)
# 
# print 'Normalized error: %g' % err2**0.5
# print 'Solvers error %g' % sol.mynorm
# print 'Congergence ratio %g' % (sol.residuals[-1]/sol.residuals[0])**(1.0/len(sol.residuals))

