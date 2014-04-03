from pylab import *
import bump_channel
import mycolors
import os

base_dir = os.path.join(os.environ['D'],'DATASTORE.RPA','projects','drag_strat')

b = bump_channel.BumpChannel(base_dir, 'taux2000_rb0110_bumplong')

lc = b.get_layers_computer()
uh = m.rdmds('')