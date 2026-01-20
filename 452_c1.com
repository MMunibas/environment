* comment
*

bomlev -1

set name 452_c1
set num 0

! Topology
OPEN UNIT 1 CARD READ NAME top_all22_prot_c1.inp
READ RTF CARD UNIT 1
CLOSE UNIT 1

OPEN UNIT 1 FORMatted READ NAME par_all22_prot_c1.inp
READ PARAmeter CARD UNIT 1
CLOSE UNIT 1

OPEN UNIT 1 CARD READ NAME @name.psf
READ PSF CARD UNIT 1
CLOSE UNIT 1

OPEN UNIT 1 CARD READ NAME @name.cor
READ COOR CARD UNIT 1
CLOSE UNIT 1

!############# ARMD ####################################################
UPDATE
OPEN UNIT 14 WRITE   FORMATTED NAME tmp_@name_@num.pdb
OPEN UNIT 9 READ FORMATTED NAME gapo_cross-c1.dat
MRMD UPAR 9 UCRG 14 PRCA 1 PRDY 1000
!#######################################################################

MMFP
GEO  sphere quartic -
     xref 0.0 yref 0.0 zref 0.0 -
     force 0.2 droff 50.0 p1 2.25 sele bynu 6:1361 end
END

coor orient MASS sele bynu 6:1361 end

WRITE COOR CARD NAME check_@name_@num.cor

wkky sele  bynu 6:1361 end
energy


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
open write unit 11 file name /data/yinc/ir/droplet_armd_kky/@name_@num_@proc.dcd
open write unit 12 file name /data/yinc/ir/droplet_armd_kky/@name_@num_@proc.vcd
open unit 13 write form name @name_@num_@proc.res

dynamics restart VERL time 0.0001 nstep 20000000 - 
iunread 10 iunwrit 13 iuncrd 11 iunvel 12 -
nprint 1000 nsavc 10 iprfrq 500 NTRFRQ 0 -
nbonds atom cdie cutnb 14.0 ctofnb 12.0 -
inbfrq -1 ihbfrq 0 IASVEL 0 ISCALE 0 

open write unit 10 card name @name_@num_@proc.pdb
write coor unit 10 pdb
close unit 10

open write unit 10 card name @name_@num_@proc.cor
write coor unit 10 card
close unit 10

stop

