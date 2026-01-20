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

stream.charmm_script('set name c1')

stream.charmm_script('set num 0')

psf_fn = "c1.psf"
read.psf_card(psf_fn)

coor_fn = "c1.cor"
read.coor_card(coor_fn, flex=True)

# Save system pdb files
write.coor_pdb("init_c1.pdb")


# Step 2: Define PhysNet energy function
#-----------------------------------------------------------

# Prepare PhysNet input parameter
selection = pycharmm.SelectAtoms(seg_id='LIG')

# Atomic numbers
ase_pept = io.read("init_c1.pdb", format="proteindatabank")
Z = ase_pept.get_atomic_numbers()
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


# Step 3: Simulation - CHARMM, PhysNet
#-----------------------------------------------------------

simulation = """

energy
mini sd nstep 750
mini abnr nstep 100

write coor card name mini_@name_@num.cor

!Heating
set proc heat
open unit 11 write form name  @name_@num_@proc.res

dynamics start time 0.0002 nstep 20000 -
firstt 48.0 finalt 300.0 teminc 10.0 ihtfrq 500 -
ieqfrq 2000 ichecw 1 twindl -5.0 twindh +5.0 iasors 0 -
nprint 100 iprfrq 500 -
!atom cdie fshift vshift cutnb 14.0 ctofnb 12.0 -
nbonds atom cdie cutnb 14.0 ctofnb 12.0 -
inbfrq -1 ihbfrq 0 -
IUNREA -1 iunwrit 11 iuncrd -1 nsavc 100 ICHECW 1


!Equilibration
open unit 10 read form name    @name_@num_@proc.res
set proc eqb
open unit 11 write form name   @name_@num_@proc.res
open write unit 20 file name   @name_@num_@proc.dcd

dynamics restart VERL time 0.0001 nstep 1000000 -
firstt 300.0 finalt 300.0 - 
nprint 1000 iprfrq 500 NTRFRQ 0 -
!atom cdie fshift vshift cutnb 14.0 ctofnb 12.0 -
nbonds atom cdie cutnb 14.0 ctofnb 12.0 -
inbfrq -1 ihbfrq 0 -
iunread 10 iunwrit 11 iuncrd 20 nsavc 1000 ICHECW 1

open write unit 10 card name @name_@num_@proc.pdb
write coor unit 10 pdb
close unit 10

open write unit 10 card name @name_@num_@proc.cor
write coor unit 10 card
close unit 10


!Production
open unit 10  read form name @name_@num_@proc.res
set proc dyna
open write unit 11 file name @name_@num_@proc.dcd
open unit 13 write form name @name_@num_@proc.res

dynamics restart VERL time 0.0001 nstep 20000000 - 
iunread 10 iunwrit 13 iuncrd 11 -
nprint 1000 nsavc 10 iprfrq 500 NTRFRQ 0 -
nbonds atom cdie cutnb 14.0 ctofnb 12.0 -
inbfrq -1 ihbfrq 0 IASVEL 0 ISCALE 0 

open write unit 10 card name @name_@num_@proc.pdb
write coor unit 10 pdb
close unit 10

open write unit 10 card name @name_@num_@proc.cor
write coor unit 10 card
close unit 10

"""
stream.charmm_script(simulation)
