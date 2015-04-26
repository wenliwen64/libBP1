gawk '{print "r "$1 $2"_a*.sac "$1 $3"_a*.sac\nrmean\nrtrend\nrotate to GCP\nw "$1"r.sac "$1"t.sac"} END {print "q"}' filelist1 | sac
dif % differentiate
int % integrate
