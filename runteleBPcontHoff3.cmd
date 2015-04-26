#!/bin/csh -f
#  runteleBPcontHoff3.cmd
#
#  UGE job for runteleBPcontHoff3 built Fri Apr 11 17:13:30 PDT 2014
#
#  The following items pertain to this script
#  Use current working directory
#$ -cwd
#  input           = /dev/null
#  output          = /u/home/l/lsmeng/cont/libBP1/runteleBPcontHoff3.joblog.$JOB_ID
#$ -o /u/home/l/lsmeng/cont/libBP1/runteleBPcontHoff3.joblog.$JOB_ID
#  error           = Merged with joblog
#$ -j y
#  The following items pertain to the user program
#  user program    = /u/home/l/lsmeng/cont/libBP1/runteleBPcontHoff3
#  arguments       = 
#  program input   = Specified by user program
#  program output  = Specified by user program
#  Resources requested
#$ -pe shared 1
#$ -l exclusive,h_data=2048M,h_rt=1:59:59
#
#  Name of application for log
#$ -v QQAPP=job
#  Email address to notify
#$ -M lsmeng@mail
#$ -m bea
#  Job is not rerunable
#$ -r n
#
# Initialization for serial execution
#
  unalias *
  set qqversion = 
  set qqapp     = "job serial"
  set qqmtasks  = 8
  set qqidir    = /u/home/l/lsmeng/cont/libBP1
  set qqjob     = runteleBPcontHoff3
  set qqodir    = /u/home/l/lsmeng/cont/libBP1
  cd     /u/home/l/lsmeng/cont/libBP1
  source /u/local/bin/qq.sge/qr.runtime
  if ($status != 0) exit (1)
#
  echo "UGE job for runteleBPcontHoff3 built Fri Apr 11 17:13:30 PDT 2014"
  echo ""
  echo "  runteleBPcontHoff3 directory:"
  echo "    "/u/home/l/lsmeng/cont/libBP1
  echo "  Submitted to UGE:"
  echo "    "$qqsubmit
  echo "  SCRATCH directory:"
  echo "    "$qqscratch
#
  echo ""
  echo "runteleBPcontHoff3 started on:   "` hostname -s `
  echo "runteleBPcontHoff3 started at:   "` date `
  echo ""
#
  source /u/local/Modules/default/init/modules.csh
  module load matlab
  setenv MCR_CACHE_ROOT $TMPDIR
#
# Run the user program
#
  echo runteleBPcontHoff3 "" \>\& runteleBPcontHoff3.output.$JOB_ID
  echo ""
  time /u/home/l/lsmeng/cont/libBP1/runteleBPcontHoff3  >& /u/home/l/lsmeng/cont/libBP1/runteleBPcontHoff3.output.$JOB_ID
#
  echo ""
  echo "runteleBPcontHoff3 finished at:  "` date `
#
# Cleanup after serial execution
#
  source /u/local/bin/qq.sge/qr.runtime
#
  echo "-------- /u/home/l/lsmeng/cont/libBP1/runteleBPcontHoff3.joblog.$JOB_ID --------" >> /u/local/apps/queue.logs/job.log.serial
 if (`wc -l /u/home/l/lsmeng/cont/libBP1/runteleBPcontHoff3.joblog.$JOB_ID  | awk '{print $1}'` >= 1000) then
        head -50 /u/home/l/lsmeng/cont/libBP1/runteleBPcontHoff3.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
        echo " "  >> /u/local/apps/queue.logs/job.log.serial
        tail -10 /u/home/l/lsmeng/cont/libBP1/runteleBPcontHoff3.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
  else
        cat /u/home/l/lsmeng/cont/libBP1/runteleBPcontHoff3.joblog.$JOB_ID >> /u/local/apps/queue.logs/job.log.serial
  endif
  exit (0)
