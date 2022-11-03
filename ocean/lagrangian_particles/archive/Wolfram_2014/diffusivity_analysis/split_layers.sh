#!/bin/bash

# functions
append_num_lines() { #{{{
  if (( $# != 2 )); then
    echo "  Usage: append_num_lines <file_name> <split num>"
  else
    echo $2 > $1_tmp.txt
    cat $1 >> $1_tmp.txt
    mv $1_tmp.txt $1
  fi
} #}}}

split_file() { #{{{
  if (( $# != 4 )); then
    echo "  Usage: split_file <in_name> <num_splits> <split_name> <append_lines_bool>"
  else
    # do the work
    numlines=`head -n 1 $1`
    #numlines=`wc $1 | awk '{print $1}'`
    #echo 'numlines = ' $numlines
    #echo 'numsplits = ' $2
    tail -n+2 $1 > $1_noheader
    # get ceil of division
    splitnum=$((($numlines + $2 - 1)/$2))
    split -a 6 -l $splitnum -d $1_noheader $1_$3
    rm $1_noheader
    if [ $4 == 'true' ] ; then
      for i in $1_$3*; do 
        append_num_lines $i $splitnum
      done
    fi
  fi
} #}}}

combine_file() { #{{{
  if (( $# != 2 )); then
    echo "  Usage: combine_file <input_name> <final_output_name>"
  else 
    # get the total number of lines and append as a header (assuming no header on files)
    #wc $1_split* | tail -n 1 | awk '{print $1}' > $2
    # assuming we have the header on the files
    head -n 1 -q $1_split* | awk '{ sum += $1 } END { print sum }' > $2
    # concatentate all the files
    for i in $1_split*; do
      tail -n+2 $i >> $2
    done
    echo 'Wrote ' $2
  fi 
} #}}}

# main 
if (( $# != 2 )); then
  echo "  Usage: split_layers.sh <cluster_data_file> <number of layers> "
  echo "  ASSUMES CURRENT DIRECTORY IS ROOT FOR SIMULATION!!!"
else 
  # split up each case in preparation for computation of gaussian kernel
  cluster_in=$1
  nlayers=$2

  DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

  # split interpolation file 
  split_file $cluster_in $nlayers layer true

fi 
#}}}
