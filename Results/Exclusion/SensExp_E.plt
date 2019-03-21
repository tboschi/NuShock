#!/data/tboschi/gnuplot/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 5.2 patchlevel 4    last modified 2018-06-01 
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2018
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
# set terminal qt 0 font "Sans,9"
# set output
unset clip points
set clip one
unset clip two
set errorbars front 1.000000 
set border 31 front lt black linewidth 1.000 dashtype solid
set zdata 
set ydata 
set xdata 
set y2data 
set x2data 
set boxwidth
set style fill   solid 0.30 noborder
set style rectangle back fc  bgnd fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02 
set style ellipse size graph 0.05, 0.03 angle 0 units xy
set dummy x, y
set format x "% h" 
set format y "e%T" 
set format x2 "% h" 
set format y2 "% h" 
set format z "% h" 
set format cb "% h" 
set format r "% h" 
set ttics format "% h"
set timefmt "%d/%m/%y,%H:%M"
set angles radians
set tics front
unset grid
unset raxis
set theta counterclockwise right
set style parallel front  lt black linewidth 2.000 dashtype solid
set key title "" center
set key at 0.0200000, 0.500000 left top vertical Left reverse enhanced autotitle nobox
set key invert samplen 1.2 spacing 1.2 width 0 height 0 
set key maxcolumns 0 maxrows 0
set key noopaque
unset label
set label 1 "color{ps191}PS191" at 0.0200000, 0.000100000, 0.00000 left norotate front textcolor rgb "#009e73"  nopoint
set label 2 "color{charm}CHARM" at 0.150000, 4.40000e-05, 0.00000 left norotate front textcolor rgb "#e51e10"  nopoint
set label 3 "color{peak}Peak search" at 0.0500000, 0.0100000, 0.00000 left norotate front textcolor rgb "#e69f00"  nopoint
set label 4 "" at 0.0200000, 0.000100000, 0.00000 left norotate back nopoint
set label 5 "Weyl state" at 0.0200000, 1.00000e-11, 0.00000 left norotate front textcolor lt -1 nopoint
set label 6 "\\shortstack{Pseudo-Dirac\\\\pair}" at 1.00000, 5.00000e-08, 0.00000 center norotate front textcolor lt -1 nopoint
set label 8 "Type I" at 0.0500000, 1.00000e-09, 0.00000 left rotate by -5 front nopoint
unset arrow
set style increment default
unset style line
unset style arrow
set style histogram clustered gap 2 title textcolor lt -1
unset object
set style textbox transparent margins  1.0,  1.0 border  lt -1 linewidth  1.0
set offsets 0, 0, 0, 0
set pointsize 1
set pointintervalbox 1
set encoding default
unset polar
unset parametric
unset decimalsign
unset micro
unset minussign
set view 60, 30, 1, 1
set view azimuth 0
set rgbmax 255
set samples 100, 100
set isosamples 10, 10
set surface 
unset contour
set cntrlabel  format '%8.3g' font '' start 5 interval 20
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5 unsorted
set cntrparam firstlinetype 0
set cntrparam points 5
set size ratio 0 1,1
set origin 0,0
set style data points
set style function lines
unset xzeroaxis
unset yzeroaxis
unset zzeroaxis
unset x2zeroaxis
unset y2zeroaxis
set xyplane relative 0.5
set tics scale  1, 0.5, 1, 1, 1
set mxtics 10.000000
set mytics 10.000000
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set mrtics default
set nomttics
set xtics border in scale 1,0.5 mirror norotate  autojustify
set xtics  norangelimit logscale autofreq 
set xtics add  (0.0500000, 0.500000, 2.00000)
set ytics border in scale 1,0.5 mirror norotate  autojustify
set ytics  norangelimit logscale 1.00000e-14,100,1.00000
set ztics border in scale 1,0.5 nomirror norotate  autojustify
set ztics  norangelimit logscale autofreq 
unset x2tics
unset y2tics
set cbtics border in scale 1,0.5 mirror norotate  autojustify
set cbtics  norangelimit logscale autofreq 
set rtics axis in scale 1,0.5 nomirror norotate  autojustify
set rtics  norangelimit autofreq 
unset ttics
set paxis 1 tics border in scale 1,0.5 nomirror norotate  autojustify
set paxis 1 tics  rangelimit autofreq 
set paxis 2 tics border in scale 1,0.5 nomirror norotate  autojustify
set paxis 2 tics  rangelimit autofreq 
set paxis 3 tics border in scale 1,0.5 nomirror norotate  autojustify
set paxis 3 tics  rangelimit autofreq 
set paxis 4 tics border in scale 1,0.5 nomirror norotate  autojustify
set paxis 4 tics  rangelimit autofreq 
set paxis 5 tics border in scale 1,0.5 nomirror norotate  autojustify
set paxis 5 tics  rangelimit autofreq 
set paxis 6 tics border in scale 1,0.5 nomirror norotate  autojustify
set paxis 6 tics  rangelimit autofreq 
set paxis 7 tics border in scale 1,0.5 nomirror norotate  autojustify
set paxis 7 tics  rangelimit autofreq 
set title "" 
set title  font "" norotate
set timestamp bottom 
set timestamp "" 
set timestamp  font "" norotate
set trange [ * : * ] noreverse nowriteback
set urange [ * : * ] noreverse nowriteback
set vrange [ * : * ] noreverse nowriteback
set xlabel "Mass (GeV)" 
set xlabel  font "" textcolor lt -1 norotate
set x2label "" 
set x2label  font "" textcolor lt -1 norotate
set xrange [ 0.0100000 : 2.00000 ] noreverse nowriteback
set x2range [ * : * ] noreverse writeback
set ylabel "$|U_{e 4}|^2$" 
set ylabel  font "" textcolor lt -1 rotate
set y2label "" 
set y2label  font "" textcolor lt -1 rotate
set yrange [ 1.00000e-12 : 1.00000 ] noreverse nowriteback
set y2range [ * : * ] noreverse writeback
set zlabel "" 
set zlabel  font "" textcolor lt -1 norotate
set zrange [ * : * ] noreverse writeback
set cblabel "" 
set cblabel  font "" textcolor lt -1 rotate
set cbrange [ * : * ] noreverse writeback
set rlabel "" 
set rlabel  font "" textcolor lt -1 norotate
set rrange [ * : * ] noreverse writeback
set paxis 1 range [ * : * ] noreverse nowriteback
set paxis 2 range [ * : * ] noreverse nowriteback
set paxis 3 range [ * : * ] noreverse nowriteback
set paxis 4 range [ * : * ] noreverse nowriteback
set paxis 5 range [ * : * ] noreverse nowriteback
set paxis 6 range [ * : * ] noreverse nowriteback
set paxis 7 range [ * : * ] noreverse nowriteback
unset logscale
set logscale z 10
set logscale y 10
set logscale x 10
set logscale cb 10
set logscale y2 10
set logscale x2 10
unset jitter
set zero 1e-08
set lmargin  -1
set bmargin  2
set rmargin  2
set tmargin  0
set locale "en_GB.UTF-8"
set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles noborder corners2color mean
set pm3d nolighting
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin screen 0.9, 0.2 size screen 0.05, 0.6 front  noinvert bdefault
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto unsorted
set loadpath 
set fontpath 
set psdir
set fit brief errorvariables nocovariancevariables errorscaling prescale nowrap v5
GNUTERM = "qt"
x = 0.0
## Last datafile plotted: "d/ExpALL_E_d_U.dat"
plot "../Limits/linetable_Weyl_E.dat"	w filledcurve lc 5 not, \
      "../Limits/linetable_TwoTwo_E.dat"	w filledcurve lc 6 not, \
      "../Limits/fadetypeI.dat"	w l lw 2 lt rgb "#ff8080" not, \
      "../Limits/linetable_Weyl_E.dat"       w l lw 2 lc 5 dt 3 not, \
      "../Limits/fulllinetable_TwoTwo_E.dat" w l lw 2 lc 6 dt 3 not, \
      "../Limits/lim_peakUe.dat"  w l not dt 1 lc 4 lw 2, \
      "../Limits/lim_ps191Ue.dat" w l not dt 1 lc 2 lw 2, \
      "../Limits/lim_charmUe.dat" w l not dt 1 lc 7 lw 2, \
      "../Limits/lim_na62Ue.dat"	   w l t "NA62" lc 8 lw 2 dt 3, \
      "../Limits/lim_shipUe.dat"	   w l t "SHiP" lc 8 lw 2 dt 2, \
      "../Limits/line_expall_e.dat" w l t "SBND" lc 8 lw 2, \
      "m/ExpALL_E_m_U.dat" w l lw 2 lc 1 t "DUNE", \
      "d/ExpALL_E_d_U.dat" w l lw 2 lc 1 dt 2 not
#    EOF
