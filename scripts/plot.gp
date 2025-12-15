
# Plot the collision integral       #
#####################################

set log x

set hidden3d

splot "OUTPUT/result.dat" using 1:2:($4 * $1**2) w l ls 2,\
      "../QCD/OUTPUT/CgXg_gg_0.txt" using 1:2:($3 / $5 * $1**2) w p ls 5 ps 0.5,\
#


splot "OUTPUT/result.dat" using 1:2:($4) w l ls 2,\
      "../QCD/OUTPUT/CgXg_gg_0.txt" using 1:2:($3 / $5) w p ls 5 ps 0.5,\
#

splot "OUTPUT/result.dat" using 1:2:($3) w l ls 3,\
      "../QCD/OUTPUT/nG0.txt" using 1:2:($3 / $4 * (2. * pi)**3) w p ls 4,\
#
