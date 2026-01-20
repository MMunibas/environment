# Test Script to import PhysNet as energy function in CHARMM via PyCHARMM

# Basics
import os
import sys
import ctypes
import pandas
import numpy as np

# ASE
from ase import Atoms
from ase import io
import ase.units as units

# PyCHARMM
import pycharmm
import pycharmm.generate as gen
import pycharmm.ic as ic
import pycharmm.coor as coor
import pycharmm.energy as energy
import pycharmm.dynamics as dyn
import pycharmm.nbonds as nbonds
import pycharmm.minimize as minimize
import pycharmm.crystal as crystal
import pycharmm.image as image
import pycharmm.psf as psf
import pycharmm.read as read
import pycharmm.write as write
import pycharmm.settings as settings
import pycharmm.lingo as stream
import pycharmm.select as select
import pycharmm.shake as shake
import pycharmm.cons_fix as cons_fix
import pycharmm.cons_harm as cons_harm
from pycharmm.lib import charmm as libcharmm
import pycharmm.lib as lib

# Step 0: Load parameter files
#-----------------------------------------------------------

# Residue and classical force field parameter

rtf_fn = "top_all22_prot_c1.inp"
read.rtf(rtf_fn)

prm_fn = "par_all22_prot_c1.inp"
read.prm(prm_fn)

settings.set_bomb_level(-2)
settings.set_warn_level(-1)

# Step 1: Read System
#-----------------------------------------------------------

stream.charmm_script('set name bulk_c1')
name = 'bulk_c1'

stream.charmm_script('set num 1')
num = '1'

psf_fn = "bulk_c1.psf"
read.psf_card(psf_fn)

coor_fn = "bulk_c1.cor"
read.coor_card(coor_fn, flex=True)

# Save system pdb files
write.coor_pdb("init_bulk_c1.pdb")

# Step 2: Set CHARMM Properties
#-----------------------------------------------------------

# Non-bonding parameter
dict_nbonds = {
    'atom': True,
    'vdw': True,
    'vswitch': True,
    'cutnb': 14,
    'ctofnb': 12,
    'ctonnb': 10,
    'cutim': 14,
    'lrc': True,
    'inbfrq': -1,
    'imgfrq': -1
    }

nbond = pycharmm.NonBondedScript(**dict_nbonds)
nbond.run()

# Custom energy
energy.show()

# PBC box
stats = coor.stat()
size = ( stats['xmax'] - stats['xmin']\
       + stats['ymax'] - stats['ymin']\
       + stats['zmax'] - stats['zmin'] ) / 3

crystal.define_cubic(length=(size + 1.0))
crystal.build(cutoff=14.0)

stream.charmm_script('image byres xcen 0.0 ycen 0.0 zcen 0.0 sele all end')

#constraints
#shake.on(bonh=True, tol=1e-7)
stream.charmm_script('wkky sele  bynu 6:6005 end')
#stream.charmm_script('cons hmcm force 1.0 refx 0.0 refy 0.0 refz 0.0 sele .not. ((segid WAT)) end')

# Deactivate CMAP correction for peptides
stream.charmm_script("skip cmap")

# Energy
energy.show()


# Step 3: Define PhysNet energy function
#-----------------------------------------------------------

# Prepare PhysNet input parameter
selection = pycharmm.SelectAtoms(seg_id='LIG')

# Atomic numbers
ase_pept = io.read("init_bulk_c1.pdb", format="proteindatabank")
Z = ase_pept.get_atomic_numbers()[:5]
print(Z)

# Checkpoint files
checkpoint = "../H2COO_12939_model_1/best_model.ckpt-4977000"

# PhysNet config file
config = "../H2COO_12939_model_1/config"

# Model units are eV and Angstrom
econv = 1./(units.kcal/units.mol)
fconv = 1./(units.kcal/units.mol)

charge = 0

# Initialize PhysNet calculator
pycharmm.MLpot(
    selection,
    fq=True,
    Z=Z,
    checkpoint=checkpoint,
    config=config,
    charge=charge,
    econv=econv,
    fconv=fconv,
    v1=True)	# Model is trained by PhysNet using tensorflow 1.x

# Custom energy
energy.show()
stream.charmm_script('cons fix sele segid LIG end')
stream.charmm_script('cons fix sele none end')
energy.show()


# Step 5: MINI - CHARMM, PhysNet
#-----------------------------------------------------------

# Do minimization or read result from last run
if True:

    stream.charmm_script('cons fix sele segid LIG end')

    # Optimization with PhysNet parameter
    minimize.run_sd(**{
        'nstep': 1000,
        'tolenr': 1e-5,
        'tolgrd': 1e-5})

    stream.charmm_script('cons fix sele none end')

    # Optimization with PhysNet parameter
    dcd_file = pycharmm.CharmmFile(
        file_name="mini_c1.dcd", file_unit=1, formatted=False, read_only=False)
    minimize.run_sd(**{
        'nstep': 1000,
        'tolenr': 1e-5,
        'tolgrd': 1e-5,
        'iuncrd': 1,
        'nsavc':  1,
        'step': 0.001})
    dcd_file.close()
    
    # Write pdb file
    write.coor_pdb("mini_c1.pdb", title="PhysNet Tripeptide in water - Minimized")
    write.coor_card("mini_c1.cor", title="PhysNet Tripeptide in water - Minimized")

else:

    # Read minimized coordinates and check energy
    read.coor_card("mini_c1.cor")
    energy.show()


# Step 4: Heating - CHARMM, PhysNet
#-----------------------------------------------------------

if True:
    
    stream.charmm_script('cons hmcm force 1.0 refx 0.0 refy 0.0 refz 0.0 sele .not. ((segid WAT)) end')
    
    timestep = 0.001	# 0.5 fs
    tottime =  5.0      # 10 ps
    savetime = 0.10     # 100 fs


    res_file = pycharmm.CharmmFile(
        file_name="heat_c1.res", file_unit=2, formatted=True, read_only=False)
    dcd_file = pycharmm.CharmmFile(
        file_name="heat_c1.dcd", file_unit=1, formatted=False, read_only=False)

    # Run some dynamics
    dynamics_dict = {
        'leap': False,
        'verlet': True,
        'cpt': False,
        'new': False,
        'langevin': False,
        'timestep': timestep,
        'start': True,
        'nstep': 2000,
        'nsavc': 100,
        'nsavv': 0,
        'inbfrq':-1,
        'ihbfrq': 50,
        'ilbfrq': 50,
        'imgfrq': 50,
        'ixtfrq': 1000,
        'iunrea':-1,
        'iunwri': res_file.file_unit,
        'iuncrd': dcd_file.file_unit,
        'nsavl':  0,  # frequency for saving lambda values in lamda-dynamics
        'iunldm':-1,
        'ilap': -1,
        'ilaf': -1,
        'nprint': 100, # Frequency to write to output
        'iprfrq': 500, # Frequency to calculate averages
        'isvfrq': 1000, # Frequency to save restart file
        'ntrfrq': 1000,
        'ihtfrq': 0,    # 200
        'ieqfrq': 1000,
        'firstt': 300,
        'finalt': 300,
        'tbath': 300,
        'iasors': 0,
        'iasvel': 1,
        'ichecw': 0,
        'iscale': 0,  # scale velocities on a restart
        'scale': 1,  # scaling factor for velocity scaling
        'echeck':-1}

    dyn_heat = pycharmm.DynamicsScript(**dynamics_dict)
    dyn_heat.run()
    
    res_file.close()
    dcd_file.close()


# Step 5: Equilibration - CHARMM, PhysNet
#-----------------------------------------------------------

if True:
        
    timestep = 0.001	# 0.2 fs
    tottime = 5.0      # 50 ps
    savetime = 0.01     # 10 fs
    
    pmass = int(np.sum(select.get_property('mass'))/50.0)
    tmass = int(pmass*10)

    str_file = pycharmm.CharmmFile(
        file_name='heat_c1.res', file_unit=3, formatted=True, read_only=False)
    res_file = pycharmm.CharmmFile(
        file_name='eqb_c1.res', file_unit=2, formatted=True, read_only=False)
    dcd_file = pycharmm.CharmmFile(
        file_name='eqb_c1.dcd', file_unit=1, formatted=False, read_only=False)

    # Run some dynamics
    dynamics_dict = {
        'leap': True,
        'verlet': False,
        'cpt': True,
        'new': False,
        'langevin': False,
        'timestep': timestep,
        'start': False,
        'restart': True,
        'nstep': 5000,
        'nsavc': 100,
        'nsavv': 0,
        'inbfrq':-1,
        'ihbfrq': 50,
        'ilbfrq': 50,
        'imgfrq': 50,
        'ixtfrq': 1000,
        'iunrea': str_file.file_unit,
        'iunwri': res_file.file_unit,
        'iuncrd': dcd_file.file_unit,
        'nsavl':  0,  # frequency for saving lambda values in lamda-dynamics
        'iunldm':-1,
        'ilap': -1,
        'ilaf': -1,
        'nprint': 1000, # Frequency to write to output
        'iprfrq': 1000, # Frequency to calculate averages
        'isvfrq': 1000, # Frequency to save restart file
        'ntrfrq': 500,
        'ihtfrq': 200,
        'ieqfrq': 0,
        'firstt': 300,
        'finalt': 300,
        'tbath': 300,
        'pint pconst pref': 1,
        'pgamma': 5,
        'pmass': pmass,
        'hoover reft': 300,
        'tmass': tmass,
        'iasors': 0,
        'iasvel': 1,
        'ichecw': 0,
        'iscale': 0,  # scale velocities on a restart
        'scale': 1,  # scaling factor for velocity scaling
        'echeck':-1}

    dyn_equi = pycharmm.DynamicsScript(**dynamics_dict)
    dyn_equi.run()
    
    str_file.close()
    res_file.close()
    dcd_file.close()

write.coor_pdb("eqb_c1.pdb")
write.coor_card("eqb_c1.cor")


# Step 6: Production - CHARMM, PhysNet
#-----------------------------------------------------------

if True:
    
    timestep = 0.0001	# 0.2 fs

    pmass = int(np.sum(select.get_property('mass'))/50.0)
    tmass = int(pmass*10)

    for ii in range(0, 10):
           
        if ii==0:

            str_file = pycharmm.CharmmFile(
                file_name='eqb_c1.res', 
                file_unit=3, formatted=True, read_only=False)
            res_file = pycharmm.CharmmFile(
                file_name='dyna_c1_{:d}.res'.format(ii), 
                file_unit=2, formatted=True, read_only=False)
            dcd_file = pycharmm.CharmmFile(
                file_name='/data/yinc/ir/bulk_physnet_kky/dyna_c1_{:d}.dcd'.format(ii), 
                file_unit=1, formatted=False, read_only=False)
            
        else:
            
            str_file = pycharmm.CharmmFile(
                file_name='dyna_c1_{:d}.res'.format(ii - 1), 
                file_unit=3, formatted=True, read_only=False)
            res_file = pycharmm.CharmmFile(
                file_name='dyna_c1_{:d}.res'.format(ii), 
                file_unit=2, formatted=True, read_only=False)
            dcd_file = pycharmm.CharmmFile(
                file_name='/data/yinc/ir/bulk_physnet_kky/dyna_c1_{:d}.dcd'.format(ii), 
                file_unit=1, formatted=False, read_only=False)
            
        # Run some dynamics
        dynamics_dict = {
            'leap': True,
            'verlet': False,
            'cpt': True,
            'new': False,
            'langevin': False,
            'timestep': timestep,
            'start': False,
            'restart': True,
#            'nstep': 100.0*1./timestep,
            'nstep': 2000000,
            'nsavc': 10,
            'nsavv': 0,
            'inbfrq':-1,
            'ihbfrq': 50,
            'ilbfrq': 50,
            'imgfrq': 50,
            'ixtfrq': 1000,
            'iunrea': str_file.file_unit,
            'iunwri': res_file.file_unit,
            'iuncrd': dcd_file.file_unit,
            'nsavl':  0,  # frequency for saving lambda values in lamda-dynamics
            'iunldm':-1,
            'ilap': -1,
            'ilaf': -1,
            'nprint': 100, # Frequency to write to output
            'iprfrq': 500, # Frequency to calculate averages
            'isvfrq': 1000, # Frequency to save restart file
            'ntrfrq': 1000,
            'ihtfrq': 0,
            'ieqfrq': 0,
            'firstt': 300,
            'finalt': 300,
            'tbath': 300,
            'pint pconst pref': 1,
            'pgamma': 5,
            'pmass': pmass,
            'hoover reft': 300,
            'tmass': tmass,
            'iasors': 0,
            'iasvel': 1,
            'ichecw': 0,
            'iscale': 0,  # scale velocities on a restart
            'scale': 1,  # scaling factor for velocity scaling
            'echeck':-1}

        dyn_prod = pycharmm.DynamicsScript(**dynamics_dict)
        dyn_prod.run()
        
        str_file.close()
        res_file.close()
        dcd_file.close()

write.coor_pdb("dyna_c1.pdb")
write.coor_card("dyna_c1.cor")

