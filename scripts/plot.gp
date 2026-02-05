#####################################
# Reference Plot script             #
#####################################

set log x

set hidden3d

splot "OUTPUT/YM.dat" using 1:2:($4 * $1**2 / (2. * pi)**3) w l ls 2,\
      "../QCD-initial/OUTPUT/CgXgg_gg_0.txt" using 1:2:($3 / $5 * $1**2) w p ls 5 ps 0.5,\
#

splot "OUTPUT/result.dat" using 1:2:($4) w l ls 2,\
      "../QCD/OUTPUT/CgXgg_gg_0.txt" using 1:2:($3 / $5 * 1.2) w p ls 5 ps 0.5,\
#

splot "OUTPUT/result.dat" using 1:2:($3 * $1**2) w l ls 3,\
      "../QCD/OUTPUT/nG0.txt" using 1:2:($3 / $4 * (2. * pi)**3 *$1**2) w p ls 4,\
#

plot "OUTPUT/result.dat" using 1:($4 * $1**2) w l ls 2,\
      "../QCD/OUTPUT/CgXgg_gg_0.txt" using 1:($3 / $5 * $1**2) w lp ls 5,\
#

set xr[:4.5]
set auto y
n=2
plot "OUTPUT/result.dat" using 1:($4 * $1**n / (2. * pi)**3) w l ls 2,\
      "../QCD/OUTPUT/CgXgg_gg_0.txt" using 1:($3 / $5 * $1**n) w lp ls 4,\
      "../QCD-initial/OUTPUT/CgXgg_gg_0.txt" using 1:($3 / $5 * $1**n) w lp ls 5,\
#

pl    "../QCD/OUTPUT/CgXgg_gg_0.txt" using 1:($3 / $5) w lp ls 4,\
      "../QCD-initial/OUTPUT/CgXgg_gg_0.txt" using 1:($3 / $5) w lp ls 5,\
#


# Hat function
set multiplot layout 2,1 title "Hat functions"
set auto
set xr[0:3]
set sample 1000; set isosample 20
splot "OUTPUT/HatQCDgg_gg.dat" u 1:2:($3 / $9 * $1**2) w pm3d,\
      exp(-x**2-y**2) * x**2
splot "OUTPUT/nG.dat" u 1:2:($3/$4 * (2 * pi)**3 * $1**2) w pm3d,\
      exp(-x**2-y**2) * x**2
unset multiplot

# Hat function
unset log x
set multiplot layout 2,1 title "Hat functions"
set xr[0:3]
set style fill transparent solid 0.5 noborder
set hidden3d
splot "OUTPUT/CgXgg_gg.dat" u 1:2:($3/$5 * (2 * pi)**3 * $1**3) w l lt -1
splot "OUTPUT/HatQCDgg_gg.dat" u 1:2:($7 / $9 * $1**3) w l lt -1
unset multiplot


set style fill transparent solid 0.5 noborder
set hidden3d
splot "data1.dat" with pm3d, \
      "data2.dat" with lines lt -1
