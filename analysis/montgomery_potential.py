from pylab import *
import os
import bump_channel
import mycolors

rcParams['figure.figsize'] = [6, 3.5]
base_dir = os.path.join(os.environ['D'],'projects','drag_strat')
figdir = './figures_variance/'

r = 'taux2000_rb0110_bump'
b = bump_channel.BumpChannel(base_dir, r)
lc = b.get_layers_computer()

T = b.rdmds('Ttave')
P = b.rdmds('Phhytave')
Piso,h = lc.interp_to_g(P, T)

Z = -h[::-1].cumsum(axis=0)[::-1]
Z = ma.masked_array(Z, Z>-5)

buoy = b.tAlpha * b.g  * lc.G
Miso = Piso - buoy * Z

#M = P - b.tAlpha * b.g * T * b.zc[:,newaxis,newaxis]
#Miso,h = lc.interp_to_g(M, T)

gidx = [21,11,6,2]

Zlevs = arange(-3000,1,200)

rc('figure.subplot', left=0.04, right=0.87, bottom=0.1, top=0.9, wspace=0.2, hspace=0.4)
figure()
for n in arange(len(gidx)):
    subplot(2,4,n+1)
    contour(Miso[gidx[n]],11, colors='k')
    title(r'M on $\theta$ = %2.1f' % lc.layers_G[gidx[n]].squeeze())
    xticks([]); yticks([]); xlabel('x'); ylabel('y')
    subplot(2,4,n+5)
    c=contourf(Z[gidx[n]], Zlevs)
    title('Depth')
    xticks([]); yticks([]); xlabel('x'); ylabel('y')
    #clabel(c, fmt='%1.1f', fontsize=7)
cb=colorbar(c, cax=axes([0.89,0.1,0.01,0.35]))

savefig('figures/M_on_theta.pdf')
