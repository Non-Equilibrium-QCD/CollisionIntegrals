#####################################
# Reference Plot script             #
#####################################

set log x

set hidden3d

splot "OUTPUT/result.dat" using 1:2:($4 * $1**2 / (2. * pi)**3) w l ls 2,\
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
