from pylab import *
from MITgcmdata import MITgcmmodel
import re
import os

flat_runs = ['taux0125_rb0110_flat',
         'taux0250_rb0110_flat',
         'taux0500_rb0110_flat',
         'taux1000_rb0110_flat',
         'taux2000_rb0110_flat',
         'taux4000_rb0110_flat',
         'taux8000_rb0110_flat' ]

bump_runs = ['taux0125_rb0110_bump',
             'taux0250_rb0110_bump',
             'taux0500_rb0110_bump',
             'taux1000_rb0110_bump',
             'taux2000_rb0110_bump',
             'taux4000_rb0110_bump',
             'taux8000_rb0110_bump']
             
double_runs = ['taux2000_rb0110_bumplong',]

drag_runs = ['taux1000_rb0014_bump',
            'taux1000_rb0028_bump', 
            'taux1000_rb0055_bump',
            'taux1000_rb0110_bump', # duplicate with bump_runs
            'taux1000_rb0220_bump',
            'taux1000_rb0440_bump',
            'taux1000_rb0880_bump']

wind_runs = flat_runs + bump_runs 
all_runs = unique(wind_runs + drag_runs + double_runs)

class BumpChannel(MITgcmmodel.ModelInstance):
    """Subclass of MITgcmdata.ModelSetup specified for the zonal channel with bump"""

    def __init__(self, base_dir, run_name, default_iter=None,
                f0=-0.8e-4, beta=1e-11, g=9.8, tAlpha=2e-4, rho0=999.8):
        
        self.f0 = f0
        self.beta = beta
        self.g = g
        self.tAlpha = tAlpha
        self.rho0 = rho0
        self.base_dir = base_dir
        
        self._set_output_dir_and_iter(base_dir, run_name, default_iter)
        MITgcmmodel.ModelInstance.__init__(self,
                output_dir=self.output_dir, grid_dir=self.grid_dir, default_iter=self.default_iter)
        self._set_physical_params()

        #self.kbottom = argmax(self.mask,axis=0)-1
        #self.kbottom[self.kbottom==0] = self.Nz-1

    def _set_output_dir_and_iter(self, base_dir, run_name, default_iter):
        self.run_name = run_name                  
        self.output_dir = os.path.join(base_dir, 'output', self.run_name)
        # figure out correct iter
        if default_iter==None:
            default_iter = 0
            for file in os.listdir(self.output_dir):
                r = re.search('Ttave.(\d{10})\.data', file)
                if r:
                    default_iter = max(default_iter, int(r.group(1)))
            print "Guessed default_iter=%10d" % default_iter
        self.default_iter = default_iter
        self.bathy = self.run_name.split('_')[-1]
        self.is_flat = (self.bathy=='flat')
        self.grid_dir = os.path.join(base_dir, 'grid_%s' % self.bathy)
        
    def _set_physical_params(self):
        # sample run name: taux1000_rb0110
        r = re.search('taux(\d+)_rb(\d+)_(bump|flat)', self.run_name)
        self.tau0 = 1e-4 * float(r.group(1))
        self.rb = 1e-5 * float(r.group(2))
        self.wind_stress = self.tau0 * sin(pi * self.yc / self.Ly)
        self.f = self.f0 + self.beta * (self.yc - self.Ly/2)

    def reinit(self, run_name, default_iter=None):
        old_grid_dir = self.grid_dir
        self._set_output_dir_and_iter(self.base_dir, run_name, default_iter)
        # reinit completly if grid_dir changed
        if old_grid_dir != self.grid_dir:
            MITgcmmodel.ModelInstance.__init__(self,
                    output_dir=self.output_dir, grid_dir=self.grid_dir, default_iter=self.default_iter)
        self._set_physical_params()
            
    def load_default_fields(self, iter=None):
        self.T = self.rdmds('Ttave', iter)
        self.V = self.rdmds('Vveltave', iter)
        self.U = self.rdmds('Uveltave', iter)
        self.W = self.rdmds('Wveltave', iter)
        #self.P = self.rdmds('Phhytave', iter)

        self.Tbar = self.T.mean(axis=2)
        self.Vbar = self.V.mean(axis=2)
        self.Ubar = self.U.mean(axis=2)
        self.Wbar = self.W.mean(axis=2)

        self.Uwave = self.U - self.Ubar[:,:,newaxis]
        self.Vwave = self.V - self.Vbar[:,:,newaxis]
        self.Wwave = self.W - self.Wbar[:,:,newaxis]
        self.Twave = self.T - self.Tbar[:,:,newaxis]

        self.UT = self.rdmds('UTtave')
        self.VT = self.rdmds('VTtave')
        self.WT = self.rdmds('WTtave')
        self.Tdif = self.rdmds('Tdiftave')

        self.VTbar = self.VT.mean(axis=2)
        self.WTbar = self.WT.mean(axis=2)
        self.Tdifbar = self.Tdif.mean(axis=2)
        
        self.UpTp = self.UT - self.U*self.cgrid_to_ugrid(self.T)
        self.UpTpbar = self.UpTp.mean(axis=2)
        self.VpTp =  self.VT - self.V*self.cgrid_to_vgrid(self.T)
        self.VpTpbar = self.VpTp.mean(axis=2)
        self.WpTp = self.WT - self.W*self.cgrid_to_wgrid(self.T)
        self.WpTpbar = self.WpTp.mean(axis=2)

        self.Tybar = self.ddy_cgrid_centered(self.Tbar)
        self.Tzbar = self.ddz_cgrid_centered(self.Tbar)
        self.s = - self.Tybar / self.Tzbar

    
    def calc_EKE(self):
        return 0.5 * ((self.rdmds('UUtave') - self.rdmds('Uveltave')**2) 
                    + (self.rdmds('VVtave') - self.rdmds('Vveltave')**2))
    
    def get_psi_iso(self):
        try:
            vflux = self.rdmds('layers_VFlux-tave', useMask=False)
        except IOError:
            vflux = self.rdmds('layers_VH-tave', useMask=False)
        return -cumsum(vflux.mean(axis=2), axis=0) * self.Lx
    
    def get_psi_iso_z(self):
        psi_iso = self.get_psi_iso()
        try:
            hv = self.rdmds('layers_HV-tave', useMask=False)
        except IOError:
            hv = self.rdmds('layers_Hs-tave', useMask=False)
        lc = self.get_layers_computer()
        return lc.transform_g_to_z(psi_iso, hv.mean(axis=2))

    def get_psi_bar(self, V=None, zpoint='C'):
        if V is None:
            V = self.rdmds('Vveltave')
        vflux = V * self.dzf_grid * self.hFacS
        psi = cumsum(vflux.mean(axis=2), axis=0) * self.Lx
        if zpoint=='F':
            return psi
        elif zpoint=='C':
            psi = vstack([zeros(self.Ny), psi])
            return 0.5 * (psi[1:] + psi[:-1])
        


