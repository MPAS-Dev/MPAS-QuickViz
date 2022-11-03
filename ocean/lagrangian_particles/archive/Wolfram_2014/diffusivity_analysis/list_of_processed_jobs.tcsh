#!/usr/bin/env tcsh

./parse_running_jobs.py | grep $1 | awk '{print $3}' | xargs -I {} basename {} | grep -o '[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}' | sort | tee $1_queue.txt
grep 'total time' $1/analyze_restart.000*/log*out* | awk '{print $0}' | xargs -I {} basename {} | grep -o '[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}' >> $1_queue.txt 

cat $1_queue.txt | sort

