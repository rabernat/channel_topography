from pylab import *

f_0=1e-4;
W_bot=1.1e-3;
delta_E=W_bot/f_0;
h_b=1e3;
ell=75e3;
L_x=2e6;
L_y=2e6;
L=1e6;
#U2=.00616
U2=arange(0.002,0.03,0.005)
Uone=ones_like(U2);
H1=1000;
#H1=[900:300/5:1200]';
H=2.985e3;
H2=H-H1;
DeltaB=1.6e-2;
kappa=.7e3;
U1mU2=DeltaB/L_y/f_0*H1;
U1=U2+U1mU2;
beta=2e-11;
gp=DeltaB/2;
F1=f_0**2 / gp / H1;
F2=f_0**2 / gp / H2;
r=W_bot / H*5;
a=U1**2 / U2 / F1 + U2/F2;
b=r / F2;
c=beta*(U1 / U2 / F1 + 1/F2);
nu1=-b / 2 / a+1j*sqrt(4*a*c-b**2)/2 / a;
nu2=-b / 2 / a-1j*sqrt(4*a*c-b**2)/2 / a;
den2=nu1-nu2;
C=-U2*f_0 / F2 / H2 / den2 / a;
dx=0.01*L;
x=arange(-L,L,dx)
xone=ones_like(x);
h=h_b*exp(-x**2/ell**2);
hx=-2*x*h/ell**2;
h=h-mean(h);
c1=cumsum(outer(Uone,h)*exp(-outer(nu1,x)), axis=1)*dx;
c2=cumsum(outer(Uone,h)*exp(-outer(nu2,x)), axis=1)*dx;
B=outer(C,xone)*exp(outer(nu1,x))*c1-(outer(C,xone)*exp(outer(nu2,x))*c2);
bx=outer(nu1,xone)*outer(C,xone)*exp(outer(nu1,x))*c1-outer(nu2,xone)*outer(C,xone)*exp(outer(nu2,x))*(c2);
barBxsquare=sum(bx**2,axis=1)*dx/L_x;
barBh_x=sum(B*outer(Uone,hx),axis=1)*dx/L_x;
tau_0=1000*f_0*barBh_x;
h1=U2**2 / (U2**2+barBxsquare)*tau_0/1000/f_0*gp/kappa/DeltaB*L_y

close('all')
rcParams['font.size'] = 8
### Wave amplitude ###
figure(1),plot(x,real(B[2,:]));grid()
#print -depsc figure1.eps

### 2d streamlines ###
y=arange(0,L_y,dx)
yone=ones(size(y));
xone=ones(size(x));
psi=-outer(y,xone*squeeze(U1[2]))+U1[2]/U2[2]*outer(yone,squeeze(B[2,:]));
figure(2, (3.5,2.8))
contour(x/1000+1000,y/1000,psi-psi.min(),10,colors='k')
xlabel('x (km)'); ylabel('y (km)')
title(r'$-U_1 y + \psi^\dagger_1$')
tight_layout()
savefig('figures/2layer/psi.pdf')
#set(gca,'fontsize',16)
#print -depsc figure2.eps

# %U=[.002,.004,.008,.012,.024];
# %tau=[.0123,.0448,.1509,.2928,.8219];%tau keeping h fixed
# %U=[.0096,.014,.0205,.0305,.047]
# %tau=[0.0499,.1010,.2014,.4007,.8071];%tau fixed h
# %U_h=[.00616,.00912,.0133,.0203,0.0312,.0510];% U varying h;
# %tau_h=[.025,.05,.1,0.1996,0.4001,.7795];%tau varying h
# %h_1=[6,478,1003,986];
# %keff=kappa*[853,16];

### tau dependence ####
figure(3),loglog(tau_0,U2,'k');grid
xlabel(r'$\tau_0$ (N m$^{-2}$)')
ylabel(r'$U_2$ (m/s)')
savez('data/U2_vs_tau_2layer.npz', tau0=tau_0, U2=U2)
#set(gca,'fontsize',16)
#print -depsc figure3.eps

show()