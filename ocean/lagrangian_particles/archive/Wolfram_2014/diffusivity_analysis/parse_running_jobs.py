#!/usr/bin/env python

import subprocess 
import numpy as np
import datetime

def job_time_left():
  """ get time left and job number """
  out = subprocess.check_output("mshow | grep pwolfram | grep Running | awk '{print $1, $5}'", shell=True)
  out = np.array(out.split()).reshape((len(out.split())/2,2))
  lefttime = np.array([(datetime.datetime.strptime(atime, '%H:%M:%S') - datetime.datetime.strptime('00:00:00', '%H:%M:%S')) for atime in out[:,1]])
  jobnum = np.array(out[:,0], dtype='int')

  return jobnum, lefttime

def get_idle_jobs():
  """ get jobnums of idle jobs"""
  out = subprocess.check_output("mshow | grep pwolfram | grep Idle | awk '{print $1}'", shell=True).split()
  jobnum = np.array(out, dtype='int')
  return jobnum

def run_time():
  """ get running time for queue"""
  runtime = subprocess.check_output("mshow | grep pwolfram | grep Running | awk '{print $1}' | xargs -I {} checkjob {} | grep WallTime | awk '{print $4}'", shell=True).split()
  runtime = np.array([(datetime.datetime.strptime(atime, '%H:%M:%S') - datetime.datetime.strptime('00:00:00', '%H:%M:%S')) for atime in runtime])

  return runtime

def job_dir(jobnum):
  """ get the job directory from run number """
  jobdir = subprocess.check_output("checkjob "+str(jobnum)+" | grep SubmitDir | awk '{print $2}'",shell=True).strip('\n')
  return jobdir

def sim_time(jobdir):
  """ get the last computed simulation time from the jobdir """
  simtimestr = subprocess.check_output("tail -n 10 "+jobdir+"/log.0000.err | grep time | tail -n 1 | awk '{print $3}' | grep -o '[0-9][0-9]_[0-9][0-9]:[0-9][0-9]:[0-9][0-9]'", shell=True).strip('\n')
  simtimestrstart = subprocess.check_output("cat "+jobdir+"/restart_timestamp | grep -o '[0-9][0-9]_[0-9][0-9]:[0-9][0-9]:[0-9][0-9]'", shell=True).strip('\n')
  simtime = datetime.datetime.strptime(simtimestr, '%d_%H:%M:%S') - datetime.datetime.strptime(simtimestrstart, '%d_%H:%M:%S')
  return simtime

def sim_end_time(jobdir):
  """ get the simulation ending time from the jobdir """
  simendtimestr = subprocess.check_output("grep config_run_duration "+jobdir+"/namelist.ocean_forward | awk '{print $3}'", shell=True).strip('\n')
  simendtime = datetime.datetime.strptime(simendtimestr.strip("'"), '%j_%H:%M:%S') - datetime.datetime(1900, 1, 1, 0, 0) + datetime.timedelta(days=1)
  return simendtime

def resubmit_job(jobnum, jobdir, submitname='384_realization.sh', walltime='16:00:00'):
  subprocess.call("canceljob " + str(jobnum), shell=True)
  subprocess.call('rm ' + jobdir + '/log.*' ,shell=True)
  subprocess.call('rm ' + jobdir + '/output*' ,shell=True)
  subprocess.call('rm ' + jobdir + '/slurm*' ,shell=True)
  subprocess.call('rm ' + jobdir + '/stats_*' ,shell=True)
  subprocess.call('mv ' + jobdir + '/restart*_backup.nc '+jobdir+'/restart*0.nc' ,shell=True)
 
  cmd = "sed -i 's/walltime=[0-9].*[0-9]/walltime="+walltime+"/' "+jobdir+"/"+submitname
  subprocess.call(cmd, shell=True)
  cmd = "cd "+jobdir+ " ; msub " + submitname
  subprocess.call(cmd, shell=True)

  return


if __name__ == "__main__":

  jobnum = get_idle_jobs()
  for ajobnum in jobnum:
    print '%d : %s waiting to run' % (ajobnum, job_dir(ajobnum))
  idlejobs = jobnum.shape[0]

  toofarbehind = 5 # percent

  jobnum, lefttime = job_time_left()
  runtime = run_time()

  okjobs = jobnum.shape[0]
  ontrack = 0
  for ajobnum, aruntime, alefttime in zip(jobnum, runtime, lefttime):
    ajobdir = job_dir(ajobnum)
    asimtime = sim_time(ajobdir)
    simtimefull = sim_end_time(ajobdir)
    percomp = asimtime.total_seconds()/simtimefull.total_seconds()*100.
    perctime = (aruntime - alefttime).total_seconds()/aruntime.total_seconds()*100.
    working = percomp > perctime
    if working:
      ontrack += 1
    print '%d : %s %20s out of %20s completed=%.2d%% time_used=%.2d%%, working=%s' % (ajobnum, ajobdir, asimtime, simtimefull, percomp, perctime, working)
    if perctime > percomp + toofarbehind:
      restart = bool(int(raw_input('Restart the job? True:1, False:0 >>> ')))
      if restart:
        print 'job %d to be restarted' % (ajobnum)
        resubmit_job(ajobnum, ajobdir)
        okjobs -= 1
  print '%d acceptable jobs are running, %d on track, %d waiting, %d total jobs' % (okjobs, ontrack, idlejobs, okjobs + idlejobs)

