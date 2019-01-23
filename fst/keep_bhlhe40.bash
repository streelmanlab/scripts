in=$1
out=$2
awk '{ if ( $2 > 7127934 && $2 < 7650069 ) print $0 }' $in > $out