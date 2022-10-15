#!/bin/bash
#for each subject you will have to modify this script for the number of runs and
#the subject number
 
#loads the fsl program
#export FSLDIR=/usr/local/bin/fsl/
.  ${FSLDIR}/etc/fslconf/fsl.sh

BASEDIR=/Users/zaz3744/Documents/current_projects/ACNlab/BrainMAPD/MID_func_conn/fsl_ppi

#sets directories and establishes what subject is being submitted for modeling
SUBJ=$1
FSLDATADIR=${BASEDIR}/fldir/${SUBJ}
TEMPLATEDIR=${BASEDIR}/templatedir
TIMECOURSEDIR=${BASEDIR}/timecourses
SEEDDIR=${BASEDIR}/seeds
DATADIR=${BASEDIR}/data
TRIMDIR=${BASEDIR}/data


#########
 
#makes the orient file
for run in 1 2; do
 #OUTPUT=${FSLDATADIR}/run${run}_output
 #echo $OUTPUT
 #makes the fsf files from the template fsf file
 cp ${TEMPLATEDIR}/design.fsf ${TEMPLATEDIR}/${SUBJ}_run${run}.fsf
 sed -i '' -e "s/10001/${SUBJ}/g" ${TEMPLATEDIR}/${SUBJ}_run${run}.fsf
 sed -i '' -e "s/run-1/run-${run}/g" ${TEMPLATEDIR}/${SUBJ}_run${run}.fsf
 sed -i '' -e "s/Run1/Run${run}/g" ${TEMPLATEDIR}/${SUBJ}_run${run}.fsf

 #need to create appropriate data files. First, cut off first two images
 #fslroi ${DATADIR}/sub-${SUBJ}_ses-2_task-MID_run-${run}_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz ${DATADIR}/sub-${SUBJ}_ses-2_task-MID_run-${run}_space-MNI152NLin6Asym_desc-preproc_bold_cut.nii.gz 2 279
 #fslmeants -i ${DATADIR}/sub-${SUBJ}_ses-2_task-MID_run-${run}_space-MNI152NLin6Asym_desc-preproc_bold_cut.nii.gz -o ${TIMECOURSEDIR}/sub-${SUBJ}_timecourse.txt -m ${SEEDDIR}/VS_8mmsphere_Oldham_Rew.nii.gz 

 #runs the analysis using the newly created fsf file
 #feat ${TEMPLATEDIR}/${SUBJ}.fsf
done
