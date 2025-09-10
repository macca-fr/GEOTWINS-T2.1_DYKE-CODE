#!/bin/csh
#plot dike paths and final shapes
gmtset LABEL_FONT_SIZE  = 35
gmtset LABEL_OFFSET     = 0.2i
gmtset ANOT_FONT_SIZE   = 30
gmtset TICK_PEN         = 2
gmtset HEADER_FONT_SIZE = 40
gmtset HEADER_OFFSET    = 0.5i
gmtset PAPER_MEDIA      = Custom_1050x900
gmtset PAGE_ORIENTATION = portrait

# ---------------------------------------------------------------------------------------------------
# -------------------- RANGE FOR THE SIMULATION DOMAIN AND PATH FOR INPUT FILES ---------------------
# ---------------------------------------------------------------------------------------------------
set Range=-R-250/250/2250/2750
set Proj=-JX10i/-10i

set outf=crack_shape.ps
set input_path=./output/
# ---------------------------------------------------------------------------------------------------
set N_dike = 1
# ---------------------------------------------------------------------------------------------------
set iter = `awk '{a=$1} END{print a}' ${input_path}M-FLUX_01.dat`
sed -i "s/D/E/g" ${input_path}xzSxxSxzSzz_r01_i${iter}.dat

# ---------------------------------------------------------------------------------------------------
# ------------ PLOT THE BACKGROUND STRESS and MAX COMPR STRESS DIRECTIONS ------------
# ---------------------------------------------------------------------------------------------------
set RES=2.5
makecpt -Ccool -D -I -T-1.5/1.5/0.5 > palette.cpt
awk '{if ($3<10 && $3>-10) print $1*1000, $2*1000, $3}' ${input_path}xzSxxSxzSzz_r01_i${iter}.dat | surface -Gtempa.grd -I$RES $Range
grdimage tempa.grd $Proj -Cpalette.cpt -K -X1.85i -Y1.5i  >  $outf

# ---------------------------------------------------------------------------------------------------
# ----------- PLOT THE DYKE PATHS AND THE DYKE SHAPES AT LAST ITERATION OF THE SIMULATION -----------
# ---------------------------------------------------------------------------------------------------

awk '{print $1*1000, $2*1000}' ${input_path}XtZtDeE_01.dat | psxy $Proj $Range -O -K -W2p,- >> $outf
awk '{print $1*1000, $2*1000}' ${input_path}crack_shape_r01_i${iter}.dat | psxy $Proj $Range -O -K -W1p -G200/200/200 -M >> $outf
# ---------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------
# ------------------------- PRINT AXIS, LABELS, TITLE AND COLOR SCALE -------------------------------
# ---------------------------------------------------------------------------------------------------
psbasemap $Proj $Range -Ba50f100":x [m]:"/a50f100":z [m]:"WSne -O -K >> $outf
psbasemap $Proj $Range -Ba50f100/a50f100N                                          -O -K >> $outf

psscale -Cpalette.cpt -D10i/5i/10i/0.5i -O -I0.1 -B0.5":@~s@~@-xx@- [MPa]:" >> $outf
# ---------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------

rm *.grd
rm *.cpt
