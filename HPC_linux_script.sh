#!/bin/sh

userCampusID=$(echo `logname`)
userProjectCheck=$(echo `whoami`)

export SUB_DIR=`pwd`
NUM_Proc=10
Queue="qGPU"
WTIMEM=$(( 10 * 1440 ))
CAMPUS_ID=`whoami`
Project=RS10263

subj_list=$(cat $1)

for s in ${subj_list} ; do

echo "#!/bin/bash">$SUB_DIR/JobSubmit-${s}.sh
echo '#SBATCH -N 1 ' >>$SUB_DIR/JobSubmit-${s}.sh
echo '#SBATCH -n ' $NUM_Proc >>$SUB_DIR/JobSubmit-${s}.sh
echo '#SBATCH -p ' $Queue  >>$SUB_DIR/JobSubmit-${s}.sh
echo '#SBATCH -t ' $WTIMEM >>$SUB_DIR/JobSubmit-${s}.sh
echo '#SBATCH -J ' $userCampusID'-'${s} >>$SUB_DIR/JobSubmit-${s}.sh
echo '#SBATCH -e '${s}'-%J.err'>>$SUB_DIR/JobSubmit-${s}.sh
echo '#SBATCH -o '${s}'-%J.out'>>$SUB_DIR/JobSubmit-${s}.sh
echo '#SBATCH -A ' $Project >>$SUB_DIR/JobSubmit-${s}.sh
echo '#SBATCH --mem-per-cpu=3000'>>$SUB_DIR/JobSubmit-${s}.sh
echo 'mkdir /runjobs/'$Project'/$SLURM_JOBID/'>> $SUB_DIR/JobSubmit-${s}.sh
echo 'cp -R '$SUB_DIR'/'${s}'.mat   /runjobs/'$Project'/$SLURM_JOBID/'>> $SUB_DIR/JobSubmit-${s}.sh
echo 'cp '$SUB_DIR'/SubmitScript-'${s}'.sh /runjobs/'$Project'/$SLURM_JOBID/' >> $SUB_DIR/JobSubmit-${s}.sh

echo 'matlab -nodisplay -r "compute /home/rshome/RS10263/WD/'${s}'.mat"' >>$SUB_DIR/SubmitScript-${s}.sh                                                                                                                                                                              

chmod 755 $SUB_DIR/SubmitScript-${s}.sh
echo '/runjobs/'$Project'/$SLURM_JOBID/SubmitScript-'${s}'.sh'>>$SUB_DIR/JobSubmit-${s}.sh
echo 'rsync -a --ignore-existing /runjobs/'$Project'/$SLURM_JOBID/  '$SUB_DIR'/'>>$SUB_DIR/JobSubmit-${s}.sh
echo 'cd /runjobs/'$Project>>$SUB_DIR/JobSubmit-${s}.sh
echo 'rm -R $SLURM_JOBID'>>$SUB_DIR/JobSubmit-${s}.sh


sbatch JobSubmit-${s}.sh
done

