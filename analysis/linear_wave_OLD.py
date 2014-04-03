from pylab import *
from scipy.special import erf, erfc

# non-dimensional parameters
eps = 0.1
h0 = 0.1
r = 1.0
l = 0.05

# coordinates
K = logspace(-2.2,1.,5000)
z = linspace(0,1,110)
x = linspace(-2*pi,4*pi,1000)

close('all')
    
def G(k, eps, l):
    return ((1.+eps-1j*l*k**-1)*k - tanh(k))/((1.+eps-1j*l*k**-1)*k*tanh(k) - 1.)
    
def calc_Phi(k,eps,h0,r=0.1,l=0.1):
    # amplitudes
    A = -eps * h0 * (-1 + 1j*k*r - eps*k*G(k,eps,l)**-1)**-1
    B = -A * G(k, eps, l)**-1
    return A*cosh(k*z) + B*sinh(k*z), A*k*sinh(k*z) + B*k*cosh(k*z)

def calc_Phi_LS(k,eps,h0,r=0.1,l=0.1):
    # binomial ap prox, wrong?
    #A = 1j*eps*h0 * (r*k + l*eps**-1*k**-1)**-1
    A = -eps*h0 / (-1 + 1j*k*r + (1 - 1j*l*(k*eps)**-1)**-1)
    B = A * (eps*k -1j*l) **-1
    return A + B*k*z, B*k
    
def calc_Psi_LS(x,a,eps,h0,r,l):
    mu1=-l/2/eps+sqrt(l**2/4/eps**2+l/r/eps)
    mu2=-l/2/eps-sqrt(l**2/4/eps**2+l/r/eps)
    print mu1,mu2
    F=-sqrt(pi/2)*a/(mu2-mu1)/r
    am1=a**2*mu1
    am2=a**2*mu2
    A = -mu1*exp(a**2*mu1**2/2+mu1*x) * (erf((x+am1)/a/sqrt(2))-1)
    A += mu2*exp(a**2*mu2**2/2+mu2*x) * (erfc(-(x+am2)/a/sqrt(2)))
    plot(A)
    intA = -exp(a**2*mu1**2/2+mu1*x) * (-1+erf((x+am1)/a/sqrt(2)))
    intA += exp(a**2*mu2**2/2+mu2*x) * (erfc(-(x+am2)/a/sqrt(2)))
    intA += -sqrt(2/pi) * exp(-x**2/2/a**2) * a / (x+am2) * (1-a**2/(x+am2)**2)
    A = F*A
    intA = F * intA
    B = l * intA + eps*A
    B = B - mean(B)
    return A[:,newaxis] * z + B[:,newaxis]

rc('figure.subplot', left=0.1, right=0.97, bottom=0.1, top=0.9, wspace=0.2, hspace=0.4)
rc('font', size=8)
if False:
    Phi = zeros((len(K),len(z)),dtype=dtype('complex128'))
    Phi_LS = zeros((len(K),len(z)),dtype=dtype('complex128'))
    for n in arange(len(K)):
        k = K[n]
        Phi[n] = calc_Phi(k,eps,h0,r,l)
        Phi_LS[n] = calc_Phi_LSi(k,eps,h0,r,l)

    LS = False
    Phi_levs = arange(0,round(abs(Phi).max()),0.5)
    close('all') 
    figure(figsize=(6.5,5.4))
    n=0
    for P in [Phi,Phi_LS]:
        n+=1
        Phi_mag = ma.masked_array(abs(P), abs(P) < 1e-6)
        # mask phases where the magnitude is zero
        Phi_phs = ma.masked_array(
            unwrap(angle(P[:,::-1]))[:,::-1], Phi_mag < 1e-6)
   
        sp1=subplot(2,2,1+2*(n-1))
        C1=contourf(K,z,Phi_mag.T/abs(eps*h0), Phi_levs, extend='max')
        contour(K,z,Phi_mag.T/abs(eps*h0), Phi_levs, colors='k')
        #clabel(C1)
        gca().set_xscale('log')
        xlabel('k')
        ylabel('z')
        title(r'$|\phi| / \epsilon h^\ast$')
        sp2=subplot(2,2,2+2*(n-1))
        C2=contourf(K,z,Phi_phs.T/pi, arange(-1,1.001,1./16), cmap=get_cmap('bwr'), extend='both')
        contour(K,z,Phi_phs.T/pi, arange(-1,1.001,1./16), colors='k')
        #clabel(C1)
        gca().set_xscale('log')
        xlabel('k')
        ylabel('z')
        title(r'$\alpha$')
        #tight_layout()
        colorbar(C1,ax=sp1)
        cb=colorbar(C2,ax=sp2, ticks=arange(-1,1.1,0.5))
        cb.set_ticklabels([r'$-\pi$', r'$-\pi/2$','0',r'$\pi/2$', r'$\pi$'])
    gcf().text(0.5,0.92,
        r'$\epsilon=$%3.2f $\lambda^\ast=$%3.2f $r^\ast=$%3.2f'% (eps,l,r),
        ha='center')
    gcf().text(0.5,0.46,'(large scale approx.)', ha='center')
    show()
    savefig('figures/linear/phi_eps-%3.2f_lam-%3.2f_r%3.2f.pdf' % (eps,l,r))
    

# fourier components of gaussian bump
N =1000 # grid points
Lx = 2000. # km
Ld = 20 # km
dx = (Lx / N) / Ld  # in n.d. units
bumpwidth = 10 * (N/400.)
x = arange(-N/2,N/2)
X = x * dx
h = exp(-x**2 / (2*bumpwidth**2))
#h = zeros(len(x)); h[N/2]=1.
K = 2*pi*fftfreq(N) / dx
hm = fft(h)

#Psi = zeros((len(z),N),dtype('complex128'))
Phi = zeros((len(z),N),dtype('complex128'))
Phi_LS = zeros_like(Phi)
PsixPsiz = zeros((len(z),N)) # heat transport
PsixPsiz_LS = zeros_like(PsixPsiz)
for n in arange(1,N):
    k = K[n]
    Phi[:,n], Phi_z = calc_Phi(k,eps,h0*hm[n],r,l)
    Phi_LS[:,n], Phi_z_LS = calc_Phi_LS(k,eps,h0*hm[n],r,l)
    PsixPsiz[:,n] = 2*real(1j*k*Phi[:,n]*conj(Phi_z))
    PsixPsiz_LS[:,n] = 2*real(1j*k*Phi_LS[:,n]*conj(Phi_z_LS))

Psi = ifft(Phi)
Psi_LS = ifft(Phi_LS)

PsixPsiz_bot = sum(PsixPsiz[0])
PsiHx_bot = sum(real(0.5*(Psi[0,1:]+Psi[0,:-1])) * diff(h))

print PsixPsiz_bot, PsiHx_bot

# Paola's analytical soln
#Psi_PC = calc_Psi_LS(X,bumpwidth,eps,h0,r,l)

if True:

    figure(figsize=(6.5, 2.0))

    psiscale = round(abs(Psi/eps/h0).max())+1
    psi_levs = arange(-psiscale -0.25 , psiscale + 0.3,0.5)
    psi_ticks = arange(-psiscale,psiscale+0.5,0.5)
    subplot(121)
    c1=contourf(X,z,real(Psi) / eps / h0, psi_levs, extend='both', cmap=get_cmap('bwr'))
    contour(X,z,real(Psi) / eps / h0, psi_levs, colors='k')
    gca().fill_between(X, h0*h, color='k')
    title(r'$\psi/\epsilon h_0$ ($\epsilon=$%3.2f $\lambda^\ast=$%3.2f $r^\ast=$%3.2f)' % (eps,l,r))
    xlabel(r'$x / L_d$')
    ylabel(r'$z/H$')
    colorbar(c1, ticks=psi_ticks)

    subplot(122)
    c1=contourf(X,z,real(Psi_LS) / eps / h0, psi_levs, extend='both', cmap=get_cmap('bwr'))
    contour(X,z,real(Psi_LS) / eps / h0, psi_levs, colors='k')
    gca().fill_between(X, h0*h, color='k')
    #title(r'$\psi/\epsilon h_0$ ($\epsilon=$%3.2f $\lambda^\ast=$%3.2f $r^\ast=$%3.2f)' % (eps,l,r))
    title("Large-Scale Approx.")
    xlabel(r'$x / L_d$')
    ylabel(r'$z/H$')
    colorbar(c1, ticks=psi_ticks)

    show()
    savefig('figures/linear/psi_eps-%3.2f_lam-%3.2f_r%3.2f.pdf' % (eps,l,r))
    
# did not work
#Psi = zeros((len(z),N),dtype('complex128'))
#for n in arange(1,N):
#    Psi += exp( 1j * K[n] * X)[newaxis,:] * Phi[:,n][:,newaxis]

