#!/bin/sh -v
#tar cvf - * | ssh -C cfd.mae.wvu.edu "bzip2 > /home/www/html/ring/dev/src/arc/domview\"$(date +%y%m%d)\"$1.tbz"
tar cvjf ../arc/sim$(date +%y%m%d)$1.tbz *.java
