#!/bin/bash

#SBATCH -A p30954
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH --mem=64G
#SBATCH -J fsl_first_level

#for each subject you will have to modify this script for the number of runs and
#the subject number
 
#loads the fsl program
#module load fsl
#.  ${FSLDIR}/etc/fslconf/fsl.sh

BASEDIR=/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/MID_func_conn/fsl_ppi

#sets directories and establishes what subject is being submitted for modeling
SUBJ=$1
FSLDATADIR=${BASEDIR}/fldir/${SUBJ}
TEMPLATEDIR=${BASEDIR}/templatedir
TIMECOURSEDIR=${BASEDIR}/timecourses
SEEDDIR=${BASEDIR}/seeds
DATADIR=${BASEDIR}/data
TRIMDIR=${BASEDIR}/data

mkdir $FSLDATADIR
#mkdir $FSLDATADIR/run-1
#mkdir $FSLDATADIR/run-2

#chmod -R 777 $FSLDATADIR

#########
 
#makes the orient file
for run in 1 2; do

# echo "creating timecourse file and trimmed functional file"
#need to create appropriate data files. First, cut off first two images
#fslroi ${DATADIR}/sub-${SUBJ}_ses-2_task-MID_run-${run}_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz ${TRIMDIR}/sub-${SUBJ}_ses-2_task-MID_run-${run}_space-MNI152NLin6Asym_desc-preproc_bold_cut.nii.gz 2 279
 fslmeants -i ${TRIMDIR}/ssub-${SUBJ}_ses-2_task-MID_run-${run}_space-MNI152NLin6Asym_desc-preproc_bold_cut.nii.gz -o ${TIMECOURSEDIR}/sub-${SUBJ}_run-${run}_timecourse.txt -m ${SEEDDIR}/VS_8mmsphere_Oldham_Rew.nii.gz 
 chmod 777 ${TIMECOURSEDIR}/sub-${SUBJ}_run-${run}_timecourse.txt
 #runs the analysis using the newly created fsf file
# echo "Starting first level model"
 feat ${TEMPLATEDIR}/${SUBJ}_run${run}.fsf
done
