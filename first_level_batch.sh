#!/bin/bash
#for each subject you will have to modify this script for the number of runs and
#the subject number
 
#loads the fsl program
#export FSLDIR=/usr/local/bin/fsl/
#.  ${FSLDIR}/etc/fslconf/fsl.sh

BASEDIR=/Users/zacharyanderson/Documents/ACNlab/BrainMAPD

#sets directories and establishes what subject is being submitted for modeling
SUBJ=$1
FSLDATADIR=${BASEDIR}/PPI/fldir/${SUBJ}
TEMPLATEDIR=${BASEDIR}/PPI
TIMECOURSEDIR=${BASEDIR}/PPI/timecourses
DATADIR=${BASEDIR}/PPI/data
SEEDDIR=${BASEDIR}/PPI/seeds

mkdir $FSLDATADIR
mkdir $FSLDATADIR/run-1
mkdir $FSLDATADIR/run-2

#########
 
#makes the orient file
for run in 1 2; do
 OUTPUT=${FSLDATADIR}/run${run}_output
 echo $OUTPUT
 #makes the fsf files from the template fsf file
 cp ${TEMPLATEDIR}/design.fsf ${TEMPLATEDIR}/${SUBJ}.fsf
 sed -i '' -e "s/10001/${SUBJ}/g" ${TEMPLATEDIR}/${SUBJ}.fsf
 sed -i '' -e "s/run-1/run-${run}/g" ${TEMPLATEDIR}/${SUBJ}.fsf
 sed -i '' -e "s/Run1/Run${run}/g" ${TEMPLATEDIR}/${SUBJ}.fsf

 #need to create appropriate data files. First, cut off first two images
 fslroi ${DATADIR}/sub-${SUBJ}_ses-2_task-MID_run-${run}_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz ${DATADIR}/sub-${SUBJ}_ses-2_task-MID_run-${run}_space-MNI152NLin6Asym_desc-preproc_bold_cut.nii.gz 2 279
 fslmeants -i ${DATADIR}/sub-${SUBJ}_ses-2_task-MID_run-${run}_space-MNI152NLin6Asym_desc-preproc_bold_cut.nii.gz -o ${TIMECOURSEDIR}/sub-${SUBJ}_timecourse.txt -m ${SEEDDIR}/VS_8mmsphere_Oldham_Rew.nii.gz 

 #runs the analysis using the newly created fsf file
 feat ${TEMPLATEDIR}/${SUBJ}.fsf
done
