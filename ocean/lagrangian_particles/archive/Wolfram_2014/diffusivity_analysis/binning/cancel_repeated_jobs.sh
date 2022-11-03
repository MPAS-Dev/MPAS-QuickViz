#!/usr/bin/env bash
# Phillip J. Wolfram
# 03/17/2016
#set -x # echo on

submitted_jobs(){
  # queued jobs
  mshow | grep pwolfram | awk '{print $1}' | xargs -I{} checkjob {} | grep SubmitDir | grep 'low' | awk '{print $2}' > queued_jobs_tmp
  cat queued_jobs_tmp | grep -o "realization_[0-9][0-9]-[0-9][0-9]" | sort > queued_jobs
  # completed jobs
  ls -al low_*_realization_*/log.0000.out.* | awk '{print $9}' | grep -o "realization.*/" > completed_jobs_tmp
  cat completed_jobs_tmp | grep -o "realization_[0-9][0-9]-[0-9][0-9]" > completed_jobs

  wc queued_jobs | awk '{print $1}' | xargs -I{} echo {}' queued jobs'
  wc completed_jobs | awk '{print $1}' | xargs -I{} echo {}' completed jobs'
  cat queued_jobs completed_jobs | wc | awk '{print $1}' | xargs -I{} echo {}' jobs completed or in progress'
  cat queued_jobs completed_jobs | sort > total_jobs

  rm -rf queued_jobs_tmp completed_jobs_tmp
}

todo_jobs(){
  ls -d -1 low_pass_realization_*/ | grep -o 'realization_[0-9][0-9]-[0-9][0-9]' > todo_jobs

}

start_unstarted(){
  for i in `comm -23 todo_jobs total_jobs`; do
    echo 'Restarting job low_pass_'$i
    cd 'low_pass_'$i
    msub 640_realization_low_pass.sh
    cd ..
  done
}

duplicated_jobs(){
  # join together to find duplicated files
  cat queued_jobs completed_jobs | sort | uniq -d > duplicatedfiles
  cat duplicatedfiles
}

cancel_duplicated_jobs(){
  mshow | grep pwolfram | awk '{print $1}' | grep 'low' \
    | xargs -I{} sh -c "echo -n {}; checkjob {} | grep SubmitDir | awk '{print \$2}'" \
    | grep -f duplicatedfiles | grep -o  '[0-9]\{3,\}' \
    | xargs -I{} canceljob {}
}

submitted_jobs
duplicated_jobs
#cancel_duplicated_jobs
#submitted_jobs

todo_jobs
start_unstarted
