#!/usr/bin/bash

#SBATCH -A p30954
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH --mem=64G
#SBATCH -J fsl_first_level

#for each subject you will have to modify this script for the number of runs and
#the subject number
 
#loads the fsl program
module load fsl
.  ${FSLDIR}/etc/fslconf/fsl.sh

BASEDIR=/projects/b1108/projects/BrainMAPD_PPI

#sets directories and establishes what subject is being submitted for modeling
SUBJ=$1
FSLDATADIR=${BASEDIR}/fldir/${SUBJ}
TEMPLATEDIR=${BASEDIR}
TIMECOURSEDIR=${BASEDIR}/timecourses
SEEDDIR=${BASEDIR}/seeds

mkdir $FSLDATADIR
mkdir $FSLDATADIR/run-1
mkdir $FSLDATADIR/run-2

#########
 
#makes the orient file
for run in 1 2; do
 DATADIR=/projects/b1108/studies/brainmapd/data/processed/neuroimaging/mid/fmriprep
 TRIMDIR=/projects/b1108/projects/BrainMAPD_PPI/data
 OUTPUT=${FSLDATADIR}/run${run}_output
 echo $OUTPUT
 #makes the fsf files from the template fsf file
 cp ${TEMPLATEDIR}/design.fsf ${TEMPLATEDIR}/${SUBJ}.fsf

 echo "changing permissions so that sed can do its thing"

 chmod 777 ${TEMPLATEDIR}/${SUBJ}.fsf

 echo "adapting template to be subject specific"
 sed -e 's/10001/${SUBJ}/' ${TEMPLATEDIR}/design.fsf > ${TEMPLATEDIR}/${SUBJ}.fsf
 sed -e 's/run-1/run-${run}/' ${TEMPLATEDIR}/design.fsf > ${TEMPLATEDIR}/${SUBJ}.fsf
 sed -e 's/Run1/Run$run/' ${TEMPLATEDIR}/design.fsf > ${TEMPLATEDIR}/${SUBJ}.fsf
 
 echo "creating timecourse file and trimmed functional file"
 #need to create appropriate data files. First, cut off first two images
 fslroi ${DATADIR}/sub-${SUBJ}/ses-2/func/sub-${SUBJ}_ses-2_task-MID_run-${run}_space-MNI152NLin6Asym_desc-preproc_bold.nii.gz ${TRIMDIR}/sub-${SUBJ}_ses-2_task-MID_run-${run}_space-MNI152NLin6Asym_desc-preproc_bold_cut.nii.gz 2 279
 fslmeants -i ${TRIMDIR}/sub-${SUBJ}_ses-2_task-MID_run-${run}_space-MNI152NLin6Asym_desc-preproc_bold_cut.nii.gz -o ${TIMECOURSEDIR}/sub-${SUBJ}_timecourse.txt -m ${SEEDDIR}/VS_8mmsphere_Oldham_Rew.nii.gz 

 #runs the analysis using the newly created fsf file
 echo "Starting first level model"
 feat ${TEMPLATEDIR}/${SUBJ}.fsf
done
