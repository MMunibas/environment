* comment
*

bomlev -1

set name asw_c2
set num 1

read rtf card name       top_all22_prot_c2.inp
read para card name      par_all22_prot_c2.inp
read psf card name       @name.psf
read coor card name      @name_@num.cor

!############# ARMD ####################################################
UPDATE
OPEN UNIT 14 WRITE   FORMATTED NAME tmp_@name_@num.pdb
OPEN UNIT 9 READ FORMATTED NAME gapo_cross-c2.dat
MRMD UPAR 9 UCRG 14 PRCA 1 PRDY 1000
!#######################################################################

!test 2 dimensions
set 6 50.0000
set 7 31.0432
set 8 31.0432

! reads in image file for translation images
crystal define rect @6 @7 @8  90. 90. 90.
crystal build cutoff 13. noper 0
IMAGE BYRES XCEN 0.0 YCEN 0.0 ZCEN 0.0 SELE ALL END

NBONDS  ATOM  VSWItch CDIE  VDW SHIFT  -
        CUTNB 13.0  CTOFNB 12.0 CTONNB 8.0  WMIN 1.5  EPS 1.0

wkky sele  bynu 9:3008 end
energy

WRITE COOR CARD NAME check_@name_@num.cor

energy
mini sd nstep 750
mini abnr nstep 100

write coor card name mini_@name_@num.cor

!Heating
set proc heat
open unit 11 write form name  @name_@num_@proc.res

DYNA STRT VERL NSTE 20000 TIME 0.0002 inbfrq 10 IMGFRQ 10 -
   IPRFRQ 100 IHTFRQ 2000 IEQFRQ 0  IHBFRQ 0 -
   IUNREA -1 IUNWRI 11 IUNCRD -1 IUNVEL -1 -
   NPRINT 100 nsavc 1000 -
   iseed 314736 475645 334837 20984 -
   FirSTT 0.0 FINALT 50 TEMINC 5 -
   TWINDH 5.0 TWINDL -5.0 -
   IASORS 1 IASVEL 1 ICHECW 0

!Equilibration
open unit 10 read form name    @name_@num_@proc.res
set proc eqb
open unit 11 write form name   @name_@num_@proc.res
open write unit 20 file name   @name_@num_@proc.dcd

DYNA RESTRT VERL NSTE 10000 TIME 0.0002 inbfrq 10 IMGFRQ 10 -
   IPRFRQ 100 IEQFRQ 2000  IHBFRQ 0 -
   IUNREA 10 IUNWRI 11 IUNCRD 20 IUNVEL -1 -
   NPRINT 100  nsavc 1000  -
   FirSTT 50 FINALT 50 -
   TWINDH 5.0 TWINDL -5.0 -
   IASORS 1 IASVEL 1 ICHECW 0

open write unit 10 card name @name_@num_@proc.pdb
write coor unit 10 pdb
close unit 10

open write unit 10 card name @name_@num_@proc.cor
write coor unit 10 card
close unit 10


!Production
open unit 10  read form name @name_@num_@proc.res
set proc dyna
open write unit 11 file name /data/yinc/ir/asw_armd_kky/@name_@num_@proc.dcd
open write unit 12 file name /data/yinc/ir/asw_armd_kky/@name_@num_@proc.vcd
open unit 13 write form name @name_@num_@proc.res

DYNA RESTART VERL NSTE 20000000  TIME 0.0001 inbfrq 10 IMGFRQ 10 -
   IPRFRQ 10000 IEQFRQ 0  IHBFRQ 0 -
   IUNREA 10 IUNWRI 13 IUNCRD 11 IUNVEL 12 -
   NPRINT 1000 nsavc 10 nsavv 10  -
   IASORS 1 IASVEL 1 ICHECW 0

open write unit 10 card name @name_@num_@proc.pdb
write coor unit 10 pdb
close unit 10

open write unit 10 card name @name_@num_@proc.cor
write coor unit 10 card
close unit 10

stop

