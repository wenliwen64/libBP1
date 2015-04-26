gawk '{if(NR==1){print "r "$1 }else{print "merge "$1 }} END{print "w append .r"}' $1 | sac
