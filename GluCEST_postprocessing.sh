#!/bin/bash

#This script executes the 7T Terra post-processing pipeline for CEST data output by the Matlab software cest2d_TERRA_SYRP. This script assumes that the following processing has been completed for a given case:
#1. The structural processing script UNI_MNI_AnatomicalParcellation.sh has been run in order to generate FAST segmentation and anatomical atlases for the 3D mp2rage UNI Image.
#2. The CEST data has been processed via the Matlab GUI cest2d_TERRA_SYRP and output in dicom format. 

#This script executes the following: dcm to nifti conversion for the Matlab generated B0, B1 and B0B1 corrected CEST files, B0 and B1 map thresholding, registration of FAST segmentation and anatomical atlases from subject 3D UNI space to subject 2D CEST space, final thresholding and masking of the CEST map,  stats extraction, and generation of QC HTML files

#######################################################################################################
##DEFINE PATHS##
strct=/import/monstrum/SYRP/7T_Terra/GluCEST/Structural #path to structural data, this is where the output of UNI_MNI_AnatomicalParcellation.sh lives 
post=/import/monstrum/SYRP/7T_Terra/GluCEST/CEST #path where processed CEST data will be output
dicoms=/import/monstrum/SYRP/7T_Terra/GluCEST/Dicoms #path to where GUI Dicoms live
SVS=/import/monstrum/SYRP/7T_Terra/SVS/01_Images/05_ROImasks #path to where the SVS ROIs are (produced by Eric), will be used for MRS v CEST comparison
#######################################################################################################

#make log file directory
if ! [ -d $post/logs ]
then
mkdir $post/logs
fi

#######################################################################################################
## IDENTIFY CASES FOR PROCESSING ##

for i in $(ls /import/monstrum/SYRP/7T_Terra/GluCEST/Dicoms/) # script will only be executed for cases in this directory
do
case=${i##*/}
echo "CASE: $case"

#check for structural data
if [ -e $strct/$case/atlases/${case}-HarvardOxford-cort.nii.gz ] && [ -e $strct/$case/atlases/${case}-HarvardOxford-sub.nii.gz ] && [ -e $strct/$case/fast/${case}_seg.nii.gz ] 
then
echo "Structural Data exists for $case"
sleep 1.5
else
echo "Oh No! Structural Data is missing. Cannot process CEST! Run UNI_MNI_AnatomicalParcellation.sh first."
sleep 1.5
fi

#check for GluCEST GUI data
if [ -d $dicoms/$case/*WASSR_B0MAP2D ] && [ -d $dicoms/$case/*B1MAP2D ] && [ -d $dicoms/$case/*B0B1CESTMAP2D ] 
then
echo "CEST GUI Data exists for $case"
sleep 1.5
else
echo "Oh No! CEST GUI Data is missing. Cannot process CEST! Analyze this case with CEST_2d_TERRA first."
sleep 1.5
fi

#determine whether case should be processed
if ! [ -d $post/$case ] && [ -d $dicoms/$case/*WASSR_B0MAP2D ] && [ -d $dicoms/$case/*B1MAP2D ] && [ -d $dicoms/$case/*B0B1CESTMAP2D ] && [ -e $strct/$case/atlases/${case}-HarvardOxford-cort.nii.gz ] && [ -e $strct/$case/atlases/${case}-HarvardOxford-sub.nii.gz ] && [ -e $strct/$case/fast/${case}_seg.nii.gz ] 
then
logfile=$post/logs/$case.log
{
echo "Running CEST Post processing pipeline on $case"
sleep 1.5
#######################################################################################################

#######################################################################################################
## CONVERT B0, B1, and B0B1 CORRECTED CEST FROM DCM TO NII ##

mkdir $post/$case
echo -e "\n-------CONVERTING DICOMS for $case-------\n\n\n"

	for seq in B0MAP B1MAP B0B1CESTMAP 
	do
	/import/speedy/scripts/melliott/dicom2nifti.sh -u -r Y $post/$case/$case-$seq.nii $dicoms/$case/S*${seq}2D/*dcm
	#-u uses dcm2nii for conversion, -r Y rescales NIFTI using Dicom Rescale/Intercept fields
	done

#######################################################################################################

#######################################################################################################
## THRESHOLD B0 AND B1 MAPS ##

echo -e "\n-------THRESHOLDING B0 for $case-------\n"

#threshold b0 from -1 to 1 (ppm shift from glutamate ppm). Note: here we are still thresholding from -1 to 1. However, we first add 10 to the raw B0 map to make all values positive, then we threshold from (-1+10)=9 to (1+10)=11. If you threshold the raw B0 from -1 to 1 and then binarize with fslmaths, fsl no longer treats 0 values as null (because negative values exist) and thus the entire FOV will be assinged a value of 1. This is the work around.
fslmaths $post/$case/$case-B0MAP.nii -add 10 $post/$case/$case-B0MAP-pos.nii.gz # make pos
fslmaths $post/$case/$case-B0MAP-pos.nii.gz -thr 9 -uthr 11 $post/$case/$case-B0MAP-thresh.nii.gz #thresh
fslmaths $post/$case/$case-B0MAP-thresh.nii.gz -bin $post/$case/$case-b0.nii.gz #binarize

#b1 includes no negative values and thus can be directly thresholded
echo -e "\n-------THRESHOLDING B1 for $case-------\n"
fslmaths $post/$case/$case-B1MAP.nii -thr 0.3 -uthr 1.3 $post/$case/$case-B1MAP-thresh.nii.gz #thresh
fslmaths $post/$case/$case-B1MAP-thresh.nii.gz -bin $post/$case/$case-b1.nii.gz #binarize 

#######################################################################################################

#######################################################################################################
## REGISTER ANATOMICAL ATLASES TO B0B1CEST ##

echo -e "\n-------CREATING 2D ATLAS for $case CEST DATA-------\n\n\n"

#extract slice from 3D UNI. Note: extract slice takes a 2D slice out of the 3D image you input that corresponds to the 2D image you input. Does not perform any registration  
gunzip $strct/$case/$case-UNI-masked.nii.gz
extract_slice.sh $strct/$case/$case-UNI-masked.nii $post/$case/$case-B0B1CESTMAP.nii $post/$case/$case-UNIslice.nii
#extract_slice.sh <3D image to extract slice from> <2D reference image> <output>

#register extracted slice to CEST with flirt normmi. Saves the rigid transform matrix so that it can be applied to atlas slices  
flirt -in $post/$case/$case-UNIslice.nii -ref $post/$case/$case-B0B1CESTMAP.nii -out $post/$case/$case-UNI-inCEST.nii -omat $post/$case/$case-UNI2CEST-matrix.mat -2D -cost normmi

#extract slice and apply transform to 3D Anatomical Atlases
echo -e "\n--------HARVAD OXFORD ATLASES-------\n"

gunzip $strct/$case/atlases/*HarvardOxford*.nii.gz

	for reg in cort sub
	do
	extract_slice.sh -MultiLabel $strct/$case/atlases/$case-HarvardOxford-$reg.nii $post/$case/$case-B0B1CESTMAP.nii $post/$case/$case-HarvardOxford-$reg-slice.nii 
	flirt -in $post/$case/$case-HarvardOxford-$reg-slice.nii -ref $post/$case/$case-B0B1CESTMAP.nii -out $post/$case/$case-2d-HarvardOxford-$reg.nii -init $post/$case/$case-UNI2CEST-matrix.mat -applyxfm -interp nearestneighbour -2D
	done

echo -e "\n-------JHU ATLAS-------\n"
gunzip $strct/$case/atlases/$case-JHU.nii.gz
extract_slice.sh  -MultiLabel $strct/$case/atlases/$case-JHU.nii $post/$case/$case-B0B1CESTMAP.nii  $post/$case/$case-JHU-slice.nii
flirt -in $post/$case/$case-JHU-slice.nii -ref $post/$case/$case-B0B1CESTMAP.nii -out $post/$case/$case-2d-JHU.nii.gz -init $post/$case/$case-UNI2CEST-matrix.mat -applyxfm -interp nearestneighbour -2D

echo -e "\n-------BNA ATLAS-------\n"
gunzip $strct/$case/atlases/$case-BNA.nii.gz
extract_slice.sh  -MultiLabel $strct/$case/atlases/$case-BNA.nii $post/$case/$case-B0B1CESTMAP.nii  $post/$case/$case-BNA-slice.nii
flirt -in $post/$case/$case-BNA-slice.nii -ref $post/$case/$case-B0B1CESTMAP.nii -out $post/$case/$case-2d-BNA.nii.gz -init $post/$case/$case-UNI2CEST-matrix.mat -applyxfm -interp nearestneighbour -2D

#extract slice from 3D FAST
echo -e "\n-------CREATING 2D SEGMENTATION for $case-------\n"
gunzip $strct/$case/fast/${case}_seg.nii.gz 
extract_slice.sh -MultiLabel $strct/$case/fast/${case}_seg.nii $post/$case/$case-B0B1CESTMAP.nii $post/$case/$case-FASTslice.nii
flirt -in $post/$case/$case-FASTslice.nii -ref $post/$case/$case-B0B1CESTMAP.nii -out $post/$case/$case-2d-FAST.nii.gz -init $post/$case/$case-UNI2CEST-matrix.mat -applyxfm -interp nearestneighbour -2D

#######################################################################################################

#######################################################################################################
## APPLY THRESHOLDED B0, B1 and GM/WM MAP (CSF removed) to B0B1CEST ##

echo -e "\n-------THRESHOLDING CEST for $case-------\n\n\n"

#threshold with b0
fslmaths $post/$case/$case-B0B1CESTMAP.nii -mul $post/$case/$case-b0.nii.gz $post/$case/$case-CEST_b0thresh.nii.gz

#threshold with b1
fslmaths $post/$case/$case-CEST_b0thresh.nii.gz -mul $post/$case/$case-b1.nii.gz $post/$case/$case-CEST_b0b1thresh.nii.gz

#threshold with GM/WM mask (also makes GM and WM maps for stats extraction)
fslmaths $post/$case/$case-2d-FAST.nii.gz -thr 2 $post/$case/$case-tissuemap.nii.gz #GM and WM, CSF removed
#CSF=1, GM=2, WM=3 in 2d FAST labelmap
fslmaths $post/$case/$case-tissuemap.nii.gz -thr 2 -uthr 2 $post/$case/$case-GM.nii.gz #GM only
fslmaths $post/$case/$case-GM.nii.gz -bin $post/$case/$case-GM.nii.gz
fslmaths $post/$case/$case-tissuemap.nii.gz -thr 3 -uthr 3 $post/$case/$case-WM.nii.gz #WM only
fslmaths $post/$case/$case-WM.nii.gz -bin $post/$case/$case-WM.nii.gz
fslmaths $post/$case/$case-tissuemap.nii.gz -bin $post/$case/$case-tissuemap-bin.nii.gz 
#perform tissue map thresholding
fslmaths $post/$case/$case-CEST_b0b1thresh.nii.gz -mul $post/$case/$case-tissuemap-bin.nii.gz $post/$case/$case-CEST-finalthresh.nii.gz 

#######################################################################################################

#######################################################################################################
## MASK THE FINAL THRESHOLDED CEST MAP. Creates a mask starting with b1, erodes, and applys to CEST map to remove non-brain voxels on image edges ##

echo -e "\n-------MASKING CEST for $case-------\n\n\n"

fslmaths $post/$case/$case-B1MAP.nii -bin $post/$case/CEST-masktmp.nii.gz
fslmaths $post/$case/CEST-masktmp.nii.gz -ero -kernel sphere 1 $post/$case/CEST-masktmp-er1.nii.gz
fslmaths $post/$case/CEST-masktmp-er1.nii.gz -ero -kernel sphere 1 $post/$case/CEST-masktmp-er2.nii.gz
fslmaths $post/$case/CEST-masktmp-er2.nii.gz -ero -kernel sphere 1 $post/$case/$case-CEST-mask.nii.gz
fslmaths $post/$case/$case-CEST-finalthresh.nii.gz -mul $post/$case/$case-CEST-mask.nii.gz $post/$case/$case-GluCEST.nii.gz #FINAL CEST OUTPUT

#######################################################################################################

#######################################################################################################
#clean up and organize

echo "quick clean up, sweep sweep"
rm -f $post/$case/*masktmp*
rm -f $post/$case/*log*
rm -f $post/$case/$case-B0MAP-pos.nii.gz 
rm -f $post/$case/$case-B0MAP-thresh.nii.gz
rm -f $post/$case/$case/-B1MAP-thresh.nii.gz
#rm -Rf $post/$case/slices

mkdir $post/$case/orig_data
mv $post/$case/$case-B0MAP.nii $post/$case/$case-B1MAP.nii $post/$case/$case-B0B1CESTMAP.nii $post/$case/orig_data
mkdir $post/$case/slices
mv $post/$case/$case*slice.nii $post/$case/slices
mkdir $post/$case/fast
mv $post/$case/$case-2d-FAST.nii.gz $post/$case/fast
mv $post/$case/$case-tissuemap.nii.gz $post/$case/fast
mv $post/$case/$case-tissuemap-bin.nii.gz $post/$case/fast
mkdir $post/$case/atlases 
mv $post/$case/$case-WM.nii.gz $post/$case/atlases
mv $post/$case/$case-GM.nii.gz $post/$case/atlases
mv $post/$case/$case*2d*.nii.gz $post/$case/atlases

gzip $strct/$case/$case-UNI-masked.nii
gzip $strct/$case/atlases/*HarvardOxford*.nii
gzip $strct/$case/fast/${case}_seg.nii
gzip $strct/$case/atlases/$case-JHU.nii
gzip $strct/$case/atlases/$case-BNA.nii
#######################################################################################################

#######################################################################################################
## CEST MEASURE EXTRACTION ##

echo -e "\n-------EXTRACTING GLUCEST VALUES FOR $case-------\n\n\n"

mkdir $post/$case/stats
cd $post/$case

3dROIstats -mask $post/$case/atlases/$case-2d-HarvardOxford-cort.nii.gz -numROI 48 -zerofill NaN -nomeanout -nzmean -nzsigma -nzvoxels -nobriklab -1DRformat $case-GluCEST.nii.gz >> $post/$case/stats/$case-HarvardOxford-cort-stats.csv #HO cortical stats

3dROIstats -mask $post/$case/atlases/$case-2d-HarvardOxford-sub.nii.gz -numROI 21 -zerofill NaN -nomeanout -nzmean -nzsigma -nzvoxels -nobriklab -1DRformat $case-GluCEST.nii.gz >> $post/$case/stats/$case-HarvardOxford-sub-stats.csv #HO subcortical stats

	for atlas in BNA GM WM
	do
	3dROIstats -mask $post/$case/atlases/$case*$atlas.nii.gz -nomeanout -nzmean -nzsigma -nzvoxels -nobriklab -1DRformat $case-GluCEST.nii.gz >> $post/$case/stats/$case-$atlas-stats.csv
	done

if ! [ -e $post/JHU-cclabels ]
then
touch $post/JHU-cclabels
echo "3 4 5 35" >> $post/JHU-cclabels
fi

if ! [ -e $post/JHU-cerebellumlabels ]
then
touch $post/JHU-cerebellumlabels
echo "1 7 15 2 9 13" >> $post/JHU-cerebellumlabels
fi

#separate output files for JHU major WM tracts and cerebellum WM 
3dROIstats -mask $post/$case/atlases/$case-2d-JHU.nii.gz -nomeanout -nzmean -nzsigma -nzvoxels -nobriklab -1DRformat -roisel ../JHU-cclabels $case-GluCEST.nii.gz >> $post/$case/stats/$case-JHU-cc-stats.csv
3dROIstats -mask $post/$case/atlases/$case-2d-JHU.nii.gz -zerofill NaN -nomeanout -nzmean -nzsigma -nzvoxels -nobriklab -1DRformat -roisel $post/JHU-cerebellumlabels $case-GluCEST.nii.gz >> $post/$case/stats/$case-JHU-cerebellum-stats.csv

cd $post

for file in $post/$case/stats/*csv
do
sed -i 's/name/Subject/g' $file 
cut -f2-3 --complement $file >> tmp.csv
mv tmp.csv $file
done

cd $post/$case/stats

#######################################################################################################

#######################################################################################################
## FORMAT CSVS ##

echo -e "\n-------FORMATTING CSVs for $case-------\n\n\n"

#Rename gray matter stats headers
sed -i 's/\<NZMean_1\>/TotalGM_mean/g' $case-GM-stats.csv
sed -i 's/\<NZcount_1\>/TotalGM_numvoxels/g' $case-GM-stats.csv
sed -i 's/\<NZSigma_1\>/TotalGM_SD/g' $case-GM-stats.csv

#Rename white matter stats headers
sed -i 's/\<NZMean_1\>/TotalWM_mean/g' $case-WM-stats.csv
sed -i 's/\<NZcount_1\>/TotalWM_numvoxels/g' $case-WM-stats.csv
sed -i 's/\<NZSigma_1\>/TotalWM_SD/g' $case-WM-stats.csv

#Rename JHU-cc white matter stats headers
sed -i 's/\<NZMean_3\>/CCgenu_mean/g' $case-JHU-cc-stats.csv
sed -i 's/\<NZcount_3\>/CCgenu_numvoxels/g' $case-JHU-cc-stats.csv
sed -i 's/\<NZSigma_3\>/CCgenu_SD/g' $case-JHU-cc-stats.csv

sed -i 's/\<NZMean_4\>/CCbody_mean/g' $case-JHU-cc-stats.csv
sed -i 's/\<NZcount_4\>/CCbody_numvoxels/g' $case-JHU-cc-stats.csv
sed -i 's/\<NZSigma_4\>/CCbody_SD/g' $case-JHU-cc-stats.csv

sed -i 's/\<NZMean_5\>/CCsplenium_mean/g' $case-JHU-cc-stats.csv
sed -i 's/\<NZcount_5\>/CCsplenium_numvoxels/g' $case-JHU-cc-stats.csv
sed -i 's/\<NZSigma_5\>/CCsplenium_SD/g' $case-JHU-cc-stats.csv

sed -i 's/\<NZMean_35\>/R_cingulum_mean/g' $case-JHU-cc-stats.csv
sed -i 's/\<NZcount_35\>/R_cingulum_numvoxels/g' $case-JHU-cc-stats.csv
sed -i 's/\<NZSigma_35\>/R_cingulum_SD/g' $case-JHU-cc-stats.csv

#Rename JHU-cerebellum white matter stats headers
sed -i 's/\<NZMean_1\>/Mid_cerpeduncle_mean/g' $case-JHU-cerebellum-stats.csv
sed -i 's/\<NZcount_1\>/Mid_cerpeduncle_numvoxels/g' $case-JHU-cerebellum-stats.csv
sed -i 's/\<NZSigma_1\>/Mid_cerpeduncle_SD/g' $case-JHU-cerebellum-stats.csv

sed -i 's/\<NZMean_2\>/Mid_cerpeduncle_pontine_mean/g' $case-JHU-cerebellum-stats.csv
sed -i 's/\<NZcount_2\>/Mid_cerpeduncle_pontine_numvoxels/g' $case-JHU-cerebellum-stats.csv
sed -i 's/\<NZSigma_2\>/Mid_cerpeduncle_pontine_SD/g' $case-JHU-cerebellum-stats.csv

sed -i 's/\<NZMean_7\>/R_corticospinal_mean/g' $case-JHU-cerebellum-stats.csv
sed -i 's/\<NZcount_7\>/R_corticospinal_numvoxels/g' $case-JHU-cerebellum-stats.csv
sed -i 's/\<NZSigma_7\>/R_corticospinal_SD/g' $case-JHU-cerebellum-stats.csv

sed -i 's/\<NZMean_9\>/R_mediallemniscus_mean/g' $case-JHU-cerebellum-stats.csv
sed -i 's/\<NZcount_9\>/R_mediallemniscus_numvoxels/g' $case-JHU-cerebellum-stats.csv
sed -i 's/\<NZSigma_9\>/R_mediallemniscus_SD/g' $case-JHU-cerebellum-stats.csv

sed -i 's/\<NZMean_13\>/Sup_cerpeduncle_mean/g' $case-JHU-cerebellum-stats.csv
sed -i 's/\<NZcount_13\>/Sup_cerpeduncle_numvoxels/g' $case-JHU-cerebellum-stats.csv
sed -i 's/\<NZSigma_13\>/Sup_cerpeduncle_SD/g' $case-JHU-cerebellum-stats.csv

sed -i 's/\<NZMean_15\>/cerpeduncle_mean/g' $case-JHU-cerebellum-stats.csv
sed -i 's/\<NZcount_15\>/cerpeduncle_numvoxels/g' $case-JHU-cerebellum-stats.csv
sed -i 's/\<NZSigma_15\>/cerpeduncle_SD/g' $case-JHU-cerebellum-stats.csv

#Rename Harvard Oxford Cortical stats headers 
sed -i 's/\<NZMean_1\>/Frontal_Pole_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_1\>/Frontal_Pole_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_1\>/Frontal_Pole_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_2\>/Insular_Cortex_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_2\>/Insular_Cortex_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_2\>/Insular_Cortex_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_3\>/SFG_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_3\>/SFG_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_3\>/SFG_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_4\>/MFG_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_4\>/MFG_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_4\>/MFG_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_5\>/IFG_parstriangularis_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_5\>/IFG_parstriangularis_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_5\>/IFG_parstriangularis_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_6\>/IFG_parsopercularis_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_6\>/IFG_parsopercularis_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_6\>/IFG_parsopercularis_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_7\>/Precentral_Gyrus_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_7\>/Precentral_Gyrus_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_7\>/Precentral_Gyrus_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_8\>/Temporal_Pole_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_8\>/Temporal_Pole_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_8\>/Temporal_Pole_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_9\>/Superior_Temporal_Gyrus_ant_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_9\>/Superior_Temporal_Gyrus_ant_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_9\>/Superior_Temporal_Gyrus_ant_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_10\>/Superior_Temporal_Gyrus_post_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_10\>/Superior_Temporal_Gyrus_post_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_10\>/Superior_Temporal_Gyrus_post_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_11\>/Middle_Temporal_Gyrus_ant_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_11\>/Middle_Temporal_Gyrus_ant_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_11\>/Middle_Temporal_Gyrus_ant_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_12\>/Middle_Temporal_Gyrus_post_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_12\>/Middle_Temporal_Gyrus_post_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_12\>/Middle_Temporal_Gyrus_post_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_13\>/Middle_Temporal_Gyrus_temporoocc_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_13\>/Middle_Temporal_Gyrus_temporoocc_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_13\>/Middle_Temporal_Gyrus_temporoocc_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_14\>/Inferior_Temporal_Gyrus_ant_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_14\>/Inferior_Temporal_Gyrus_ant_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_14\>/Inferior_Temporal_Gyrus_ant_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_15\>/Inferior_Temporal_Gyrus_post_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_15\>/Inferior_Temporal_Gyrus_post_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_15\>/Inferior_Temporal_Gyrus_post_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_16\>/Inferior_Temporal_Gyrus_temporocc_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_16\>/Inferior_Temporal_Gyrus_temporocc_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_16\>/Inferior_Temporal_Gyrus_temporocc_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_17\>/Postcentral_Gyrus_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_17\>/Postcentral_Gyrus_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_17\>/Postcentral_Gyrus_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_18\>/Superior_Parietal_Lobule_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_18\>/Superior_Parietal_Lobule_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_18\>/Superior_Parietal_Lobule_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_19\>/Supramarginal_Gyrus_ant_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_19\>/Supramarginal_Gyrus_ant_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_19\>/Supramarginal_Gyrus_ant_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_20\>/Supramarginal_Gyrus_post_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_20\>/Supramarginal_Gyrus_post_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_20\>/Supramarginal_Gyrus_post_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_21\>/Angular_Gyrus_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_21\>/Angular_Gyrus_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_21\>/Angular_Gyrus_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_22\>/Lateral_Occipital_Cortex_sup_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_22\>/Lateral_Occipital_Cortex_sup_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_22\>/Lateral_Occipital_Cortex_sup_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_23\>/Lateral_Occipital_Cortex_inf_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_23\>/Lateral_Occipital_Cortex_inf_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_23\>/Lateral_Occipital_Cortex_inf_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_24\>/Intracalcarine_Cortex_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_24\>/Intracalcarine_Cortex_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_24\>/Intracalcarine_Cortex_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_25\>/Frontal_Medial_Cortex_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_25\>/Frontal_Medial_Cortex_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_25\>/Frontal_Medial_Cortex_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_26\>/Juxtapositional_Lobule_Cortex_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_26\>/Juxtapositional_Lobule_Cortex_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_26\>/Juxtapositional_Lobule_Cortex_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_27\>/Subcallosal_Cortex_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_27\>/Subcallosal_Cortex_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_27\>/Subcallosal_Cortex_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_28\>/Paracingulate_Gyrus_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_28\>/Paracingulate_Gyrus_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_28\>/Paracingulate_Gyrus_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_29\>/Anterior_cingulate_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_29\>/Anterior_cingulate_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_29\>/Anterior_cingulate_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_30\>/Posterior_cingulate_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_30\>/Posterior_cingulate_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_30\>/Posterior_cingulate_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_31\>/Precuneous_Cortex_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_31\>/Precuneous_Cortex_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_31\>/Precuneous_Cortex_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_32\>/Cuneal_Cortex_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_32\>/Cuneal_Cortex_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_32\>/Cuneal_Cortex_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_33\>/OFC_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_33\>/OFC_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_33\>/OFC_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_34\>/Parahippocampal_Gyrus_ant_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_34\>/Parahippocampal_Gyrus_ant_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_34\>/Parahippocampal_Gyrus_ant_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_35\>/Parahippocampal_Gyrus_post_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_35\>/Parahippocampal_Gyrus_post_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_35\>/Parahippocampal_Gyrus_post_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_36\>/Lingual_Gyrus_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_36\>/Lingual_Gyrus_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_36\>/Lingual_Gyrus_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_37\>/Temporal_Fusiform_Cortex_ant_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_37\>/Temporal_Fusiform_Cortex_ant_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_37\>/Temporal_Fusiform_Cortex_ant_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_38\>/Temporal_Fusiform_Cortex_post_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_38\>/Temporal_Fusiform_Cortex_post_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_38\>/Temporal_Fusiform_Cortex_post_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_39\>/Temporal_Occipital_Fusiform_Cortex_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_39\>/Temporal_Occipital_Fusiform_Cortex_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_39\>/Temporal_Occipital_Fusiform_Cortex_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_40\>/Occipital_Fusiform_Gyrus_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_40\>/Occipital_Fusiform_Gyrus_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_40\>/Occipital_Fusiform_Gyrus_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_41\>/Frontal_Operculum_Cortex_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_41\>/Frontal_Operculum_Cortex_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_41\>/Frontal_Operculum_Cortex_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_42\>/Central_Opercular_Cortex_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_42\>/Central_Opercular_Cortex_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_42\>/Central_Opercular_Cortex_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_43\>/Parietal_Operculum_Cortex_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_43\>/Parietal_Operculum_Cortex_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_43\>/Parietal_Operculum_Cortex_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_44\>/Planum_Polare_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_44\>/Planum_Polare_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_44\>/Planum_Polare_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_45\>/Heschls_Gyrus_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_45\>/Heschls_Gyrus_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_45\>/Heschls_Gyrus_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_46\>/Planum_Temporale_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_46\>/Planum_Temporale_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_46\>/Planum_Temporale_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_47\>/Supracalcarine_Cortex_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_47\>/Supracalcarine_Cortex_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_47\>/Supracalcarine_Cortex_SD/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZMean_48\>/Occipital_Pole_mean/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZcount_48\>/Occipital_Pole_numvoxels/g' $case-HarvardOxford-cort-stats.csv
sed -i 's/\<NZSigma_48\>/Occipital_Pole_SD/g' $case-HarvardOxford-cort-stats.csv

#Rename Harvard Oxford Subcortical stats headers 
sed -i 's/\<NZMean_1\>/Left_Cerebral_White_Matter_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_1\>/Left_Cerebral_White_Matter_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_1\>/Left_Cerebral_White_Matter_SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_2\>/Left_Cerebral_Cortex__mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_2\>/Left_Cerebral_Cortex__numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_2\>/Left_Cerebral_Cortex__SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_3\>/Left_Lateral_Ventrical_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_3\>/Left_Lateral_Ventrical_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_3\>/Left_Lateral_Ventrical_SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_4\>/Left_Thalamus_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_4\>/Left_Thalamus_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_4\>/Left_Thalamus_SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_5\>/Left_Caudate_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_5\>/Left_Caudate_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_5\>/Left_Caudate_SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_6\>/Left_Putamen_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_6\>/Left_Putamen_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_6\>/Left_Putamen_SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_7\>/Left_Pallidum_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_7\>/Left_Pallidum_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_7\>/Left_Pallidum_SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_8\>/Brain-Stem_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_8\>/Brain-Stem_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_8\>/Brain-Stem_SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_9\>/Left_Hippocampus_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_9\>/Left_Hippocampus_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_9\>/Left_Hippocampus_SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_10\>/Left_Amygdala_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_10\>/Left_Amygdala_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_10\>/Left_Amygdala_SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_11\>/Left_Accumbens_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_11\>/Left_Accumbens_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_11\>/Left_Accumbens_SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_12\>/Right_Cerebral_White_Matter_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_12\>/Right_Cerebral_White_Matter_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_12\>/Right_Cerebral_White_Matter_SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_13\>/Right_Cerebral_Cortex__mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_13\>/Right_Cerebral_Cortex__numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_13\>/Right_Cerebral_Cortex__SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_14\>/Right_Lateral_Ventricle_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_14\>/Right_Lateral_Ventricle_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_14\>/Right_Lateral_Ventricle_SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_15\>/Right_Thalamus_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_15\>/Right_Thalamus_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_15\>/Right_Thalamus_SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_16\>/Right_Caudate_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_16\>/Right_Caudate_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_16\>/Right_Caudate_SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_17\>/Right_Putamen_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_17\>/Right_Putamen_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_17\>/Right_Putamen_SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_18\>/Right_Pallidum_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_18\>/Right_Pallidum_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_18\>/Right_Pallidum_SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_19\>/Right_Hippocampus_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_19\>/Right_Hippocampus_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_19\>/Right_Hippocampus_SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_20\>/Right_Amygdala_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_20\>/Right_Amygdala_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_20\>/Right_Amygdala_SD/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZMean_21\>/Right_Accumbens_mean/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZcount_21\>/Right_Accumbens_numvoxels/g' $case-HarvardOxford-sub-stats.csv
sed -i 's/\<NZSigma_21\>/Right_Accumbens_SD/g' $case-HarvardOxford-sub-stats.csv

if ! [ -e $post/GluCEST-HarvardOxford-Cortical-Measures.csv ]
then
touch $post/GluCEST-HarvardOxford-Cortical-Measures.csv
echo "Subject	Frontal_Pole_mean	Frontal_Pole_numvoxels	Frontal_Pole_SD	Insular_Cortex_mean	Insular_Cortex_numvoxels	Insular_Cortex_SD	SFG_mean	SFG_numvoxels	SFG_SD	MFG_mean	MFG_numvoxels	MFG_SD	IFG_parstriangularis_mean	IFG_parstriangularis_numvoxels	IFG_parstriangularis_SD	IFG_parsopercularis_mean	IFG_parsopercularis_numvoxels	IFG_parsopercularis_SD	Precentral_Gyrus_mean	Precentral_Gyrus_numvoxels	Precentral_Gyrus_SD	Temporal_Pole_mean	Temporal_Pole_numvoxels	Temporal_Pole_SD	Superior_Temporal_Gyrus_ant_mean	Superior_Temporal_Gyrus_ant_numvoxels	Superior_Temporal_Gyrus_ant_SD	Superior_Temporal_Gyrus_post_mean	Superior_Temporal_Gyrus_post_numvoxels	Superior_Temporal_Gyrus_post_SD	Middle_Temporal_Gyrus_ant_mean	Middle_Temporal_Gyrus_ant_numvoxels	Middle_Temporal_Gyrus_ant_SD	Middle_Temporal_Gyrus_post_mean	Middle_Temporal_Gyrus_post_numvoxels	Middle_Temporal_Gyrus_post_SD	Middle_Temporal_Gyrus_temporoocc_mean	Middle_Temporal_Gyrus_temporoocc_numvoxels	Middle_Temporal_Gyrus_temporoocc_SD	Inferior_Temporal_Gyrus_ant_mean	Inferior_Temporal_Gyrus_ant_numvoxels	Inferior_Temporal_Gyrus_ant_SD	Inferior_Temporal_Gyrus_post_mean	Inferior_Temporal_Gyrus_post_numvoxels	Inferior_Temporal_Gyrus_post_SD	Inferior_Temporal_Gyrus_temporocc_mean	Inferior_Temporal_Gyrus_temporocc_numvoxels	Inferior_Temporal_Gyrus_temporocc_SD	Postcentral_Gyrus_mean	Postcentral_Gyrus_numvoxels	Postcentral_Gyrus_SD	Superior_Parietal_Lobule_mean	Superior_Parietal_Lobule_numvoxels	Superior_Parietal_Lobule_SD	Supramarginal_Gyrus_ant_mean	Supramarginal_Gyrus_ant_numvoxels	Supramarginal_Gyrus_ant_SD	Supramarginal_Gyrus_post_mean	Supramarginal_Gyrus_post_numvoxels	Supramarginal_Gyrus_post_SD	Angular_Gyrus_mean	Angular_Gyrus_numvoxels	Angular_Gyrus_SD	Lateral_Occipital_Cortex_sup_mean	Lateral_Occipital_Cortex_sup_numvoxels	Lateral_Occipital_Cortex_sup_SD	Lateral_Occipital_Cortex_inf_mean	Lateral_Occipital_Cortex_inf_numvoxels	Lateral_Occipital_Cortex_inf_SD	Intracalcarine_Cortex_mean	Intracalcarine_Cortex_numvoxels	Intracalcarine_Cortex_SD	Frontal_Medial_Cortex_mean	Frontal_Medial_Cortex_numvoxels	Frontal_Medial_Cortex_SD	Juxtapositional_Lobule_Cortex_mean	Juxtapositional_Lobule_Cortex_numvoxels	Juxtapositional_Lobule_Cortex_SD	Subcallosal_Cortex_mean	Subcallosal_Cortex_numvoxels	Subcallosal_Cortex_SD	Paracingulate_Gyrus_mean	Paracingulate_Gyrus_numvoxels	Paracingulate_Gyrus_SD	Anterior_cingulate_mean	Anterior_cingulate_numvoxels	Anterior_cingulate_SD	Posterior_cingulate_mean	Posterior_cingulate_numvoxels	Posterior_cingulate_SD	Precuneous_Cortex_mean	Precuneous_Cortex_numvoxels	Precuneous_Cortex_SD	Cuneal_Cortex_mean	Cuneal_Cortex_numvoxels	Cuneal_Cortex_SD	OFC_mean	OFC_numvoxels	OFC_SD	Parahippocampal_Gyrus_ant_mean	Parahippocampal_Gyrus_ant_numvoxels	Parahippocampal_Gyrus_ant_SD	Parahippocampal_Gyrus_post_mean	Parahippocampal_Gyrus_post_numvoxels	Parahippocampal_Gyrus_post_SD	Lingual_Gyrus_mean	Lingual_Gyrus_numvoxels	Lingual_Gyrus_SD	Temporal_Fusiform_Cortex_ant_mean	Temporal_Fusiform_Cortex_ant_numvoxels	Temporal_Fusiform_Cortex_ant_SD	Temporal_Fusiform_Cortex_post_mean	Temporal_Fusiform_Cortex_post_numvoxels	Temporal_Fusiform_Cortex_post_SD	Temporal_Occipital_Fusiform_Cortex_mean	Temporal_Occipital_Fusiform_Cortex_numvoxels	Temporal_Occipital_Fusiform_Cortex_SD	Occipital_Fusiform_Gyrus_mean	Occipital_Fusiform_Gyrus_numvoxels	Occipital_Fusiform_Gyrus_SD	Frontal_Operculum_Cortex_mean	Frontal_Operculum_Cortex_numvoxels	Frontal_Operculum_Cortex_SD	Central_Opercular_Cortex_mean	Central_Opercular_Cortex_numvoxels	Central_Opercular_Cortex_SD	Parietal_Operculum_Cortex_mean	Parietal_Operculum_Cortex_numvoxels	Parietal_Operculum_Cortex_SD	Planum_Polare_mean	Planum_Polare_numvoxels	Planum_Polare_SD	Heschls_Gyrus_mean	Heschls_Gyrus_numvoxels	Heschls_Gyrus_SD	Planum_Temporale_mean	Planum_Temporale_numvoxels	Planum_Temporale_SD	Supracalcarine_Cortex_mean	Supracalcarine_Cortex_numvoxels	Supracalcarine_Cortex_SD	Occipital_Pole_mean	Occipital_Pole_numvoxels	Occipital_Pole_SD" >> $post/GluCEST-HarvardOxford-Cortical-Measures.csv
fi

sed -n "2p" $post/$case/stats/$case-HarvardOxford-cort-stats.csv >> $post/GluCEST-HarvardOxford-Cortical-Measures.csv

if ! [ -e $post/GluCEST-HarvardOxford-Subcortical-Measures.csv ]
then
touch $post/GluCEST-HarvardOxford-Subcortical-Measures.csv
echo "Subject	Left_Cerebral_White_Matter_mean	Left_Cerebral_White_Matter_numvoxels	Left_Cerebral_White_Matter_SD	Left_Cerebral_Cortex__mean	Left_Cerebral_Cortex__numvoxels	Left_Cerebral_Cortex__SD	Left_Lateral_Ventrical_mean	Left_Lateral_Ventrical_numvoxels	Left_Lateral_Ventrical_SD	Left_Thalamus_mean	Left_Thalamus_numvoxels	Left_Thalamus_SD	Left_Caudate_mean	Left_Caudate_numvoxels	Left_Caudate_SD	Left_Putamen_mean	Left_Putamen_numvoxels	Left_Putamen_SD	Left_Pallidum_mean	Left_Pallidum_numvoxels	Left_Pallidum_SD	Brain-Stem_mean	Brain-Stem_numvoxels	Brain-Stem_SD	Left_Hippocampus_mean	Left_Hippocampus_numvoxels	Left_Hippocampus_SD	Left_Amygdala_mean	Left_Amygdala_numvoxels	Left_Amygdala_SD	Left_Accumbens_mean	Left_Accumbens_numvoxels	Left_Accumbens_SD	Right_Cerebral_White_Matter_mean	Right_Cerebral_White_Matter_numvoxels	Right_Cerebral_White_Matter_SD	Right_Cerebral_Cortex__mean	Right_Cerebral_Cortex__numvoxels	Right_Cerebral_Cortex__SD	Right_Lateral_Ventricle_mean	Right_Lateral_Ventricle_numvoxels	Right_Lateral_Ventricle_SD	Right_Thalamus_mean	Right_Thalamus_numvoxels	Right_Thalamus_SD	Right_Caudate_mean	Right_Caudate_numvoxels	Right_Caudate_SD	Right_Putamen_mean	Right_Putamen_numvoxels	Right_Putamen_SD	Right_Pallidum_mean	Right_Pallidum_numvoxels	Right_Pallidum_SD	Right_Hippocampus_mean	Right_Hippocampus_numvoxels	Right_Hippocampus_SD	Right_Amygdala_mean	Right_Amygdala_numvoxels	Right_Amygdala_SD	Right_Accumbens_mean	Right_Accumbens_numvoxels	Right_Accumbens_SD" >> $post/GluCEST-HarvardOxford-Subcortical-Measures.csv
fi

sed -n "2p" $post/$case/stats/$case-HarvardOxford-sub-stats.csv >> $post/GluCEST-HarvardOxford-Subcortical-Measures.csv

if ! [ -e $post/GluCEST-TotalGrayMatter-Measures.csv ]
then
touch $post/GluCEST-TotalGrayMatter-Measures.csv
echo "Subject	TotalGM_mean	TotalGM_numvoxels	TotalGM_SD" >> $post/GluCEST-TotalGrayMatter-Measures.csv
fi

sed -n "2p" $post/$case/stats/$case-GM-stats.csv >> $post/GluCEST-TotalGrayMatter-Measures.csv

if ! [ -e $post/GluCEST-TotalWhiteMatter-Measures.csv ]
then
touch $post/GluCEST-TotalWhiteMatter-Measures.csv
echo "Subject	TotalWM_mean	TotalWM_numvoxels	TotalWM_SD" >> $post/GluCEST-TotalWhiteMatter-Measures.csv
fi

sed -n "2p" $post/$case/stats/$case-WM-stats.csv >> $post/GluCEST-TotalWhiteMatter-Measures.csv

if ! [ -e $post/GluCEST-JHU-WM-Measures.csv ] 
then
touch $post/GluCEST-JHU-WM-Measures.csv
echo "Subject	CCgenu_mean	CCgenu_numvoxels	CCgenu_SD	CCbody_mean	CCbody_numvoxels	CCbody_SD	CCsplenium_mean	CCsplenium_numvoxels	CCsplenium_SD	R_cingulum_mean	R_cingulum_numvoxels	R_cingulum_SD" >> $post/GluCEST-JHU-WM-Measures.csv
fi

sed -n "2p" $post/$case/stats/$case-JHU-cc-stats.csv >> $post/GluCEST-JHU-WM-Measures.csv

if ! [ -e $post/GluCEST-Cerebellum-WM-Measures.csv ]
then
touch $post/GluCEST-Cerebellum-WM-Measures.csv
echo "Subject	Mid_cerpeduncle_mean	Mid_cerpeduncle_numvoxels	Mid_cerpeduncle_SD	Mid_cerpeduncle_pontine_mean	Mid_cerpeduncle_pontine_numvoxels	Mid_cerpeduncle_pontine_SD	R_corticospinal_mean	R_corticospinal_numvoxels	R_corticospinal_SD	R_mediallemniscus_mean	R_mediallemniscus_numvoxels	R_mediallemniscus_SD	Sup_cerpeduncle_mean	Sup_cerpeduncle_numvoxels	Sup_cerpeduncle_SD	cerpeduncle_mean	cerpeduncle_numvoxels	cerpeduncle_SD" >> $post/GluCEST-Cerebellum-WM-Measures.csv
fi

sed -n "2p" $post/$case/stats/$case-JHU-cerebellum-stats.csv >> $post/GluCEST-Cerebellum-WM-Measures.csv
#######################################################################################################

#######################################################################################################
## MAKE QC IMAGES FOR GUI OUTPUT (B0, B1, BOB1CEST) and PROCESSED CEST/ATLASES and CORRESPONDING HTMLS ##

echo -e "\n-------MAKING QC IMAGES for $case-------\n"

if ! [ -d $post/Quality_Control ]
then
mkdir $post/Quality_Control
fi

if ! [ -d $post/Quality_Control/QC_GUI ]
then
mkdir $post/Quality_Control/QC_GUI
fi

if ! [ -e $post/Quality_Control/QC_GUI/QC_GUI_CEST.html ]
then
touch $post/Quality_Control/QC_GUI/QC_GUI_CEST.html
echo "<html>

<head>
<title>B0B1 Corrected CEST QC</title>
</head>

<body>

<br><strong><font size="7.5" color="1814A1">7T TERRA CEST DATA QUALITY CONTROL</font></strong><br>
<br>
<br><strong><font size="5" color="1814A1">&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;B0 &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; B1 &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&nbsp;&nbsp;B0B1-corrected CEST</font></strong><br>
<br>

" >> $post/Quality_Control/QC_GUI/QC_GUI_CEST.html
fi

slices $post/$case/orig_data/$case-B0MAP.nii -i -1.3 1.3 -o $post/Quality_Control/QC_GUI/$case-B0MAP-qc
slices $post/$case/orig_data/$case-B1MAP.nii -i 0 3.3 -o $post/Quality_Control/QC_GUI/$case-B1MAP-qc
slices $post/$case/orig_data/$case-B0B1CESTMAP.nii -i 0 16 -o $post/Quality_Control/QC_GUI/$case-CEST-qc


echo "<br> <strong>$case</strong><br>
<img src=$case-B0MAP-qc height="350" width="350">
<img src=$case-B1MAP-qc height="350" width="350">
<img src=$case-CEST-qc height="350" width="350">
<br> 

" >> $post/Quality_Control/QC_GUI/QC_GUI_CEST.html
 

if ! [ -d $post/Quality_Control/QC_Postprocessing ]
then
mkdir $post/Quality_Control/QC_Postprocessing
fi

if ! [ -d $post/Quality_Control/QC_Postprocessing ]
then
mkdir $post/Quality_Control/QC_Postprocessing
fi

if ! [ -e $post/Quality_Control/QC_Postprocessing/QC_GluCEST.html ]
then
touch $post/Quality_Control/QC_Postprocessing/QC_GluCEST.html
echo "<html>

<head>
<title>Quality Control: FINAL GLUCEST</title>
</head>

<body>

<br><strong><font size="7.5" color="1814A1">QUALITY CONTROL: FINAL THRESHOLDED GLUCEST</font></strong><br>
<br>

" >> $post/Quality_Control/QC_Postprocessing/QC_GluCEST.html
fi

if ! [ -e $post/Quality_Control/QC_Postprocessing/QC_Atlases_Cortical.html ]
then
touch $post/Quality_Control/QC_Postprocessing/QC_Atlases_Cortical.html
echo "<html>

<head>
<title>Quality Control: GLUCEST CORTICAL ATLASES</title>
</head>

<body>

<br><strong><font size="7.5" color="1814A1">QUALITY CONTROL: GLUCEST CORTICAL ATLASES</font></strong><br>
<br>

" >> $post/Quality_Control/QC_Postprocessing/QC_Atlases_Cortical.html
fi

if ! [ -e $post/Quality_Control/QC_Postprocessing/QC_Atlases_WM.html ]
then
touch $post/Quality_Control/QC_Postprocessing/QC_Atlases_WM.html
echo "<html>

<head>
<title>Quality Control: GLUCEST WM ATLASES</title>
</head>

<body>

<br><strong><font size="7.5" color="1814A1">QUALITY CONTROL: GLUCEST WHITE MATTER ATLASES</font></strong><br>
<br>

" >> $post/Quality_Control/QC_Postprocessing/QC_Atlases_WM.html
fi

slices $post/$case/${case}-GluCEST.nii.gz -i 0 16 -o $post/Quality_Control/QC_Postprocessing/$case-GluCEST-qc
overlay 1 0 $post/$case/${case}-GluCEST.nii.gz -a $post/$case/atlases/$case-2d-HarvardOxford-cort.nii.gz 1 20 $post/Quality_Control/QC_Postprocessing/$case-HO-overlay.nii.gz 
slicer $post/Quality_Control/QC_Postprocessing/$case-HO-overlay.nii.gz  -l /import/monstrum/Applications/fsl/etc/luts/renderhsv.lut -z 0.5 $post/Quality_Control/QC_Postprocessing/$case-atlas-qc.png
rm -f $post/Quality_Control/QC_Postprocessing/$case-HO-overlay.nii.gz
overlay 1 0 $post/$case/${case}-GluCEST.nii.gz -a $post/$case/atlases/$case-2d-JHU.nii.gz 1 20 $post/Quality_Control/QC_Postprocessing/$case-JHU-overlay.nii.gz
slicer $post/Quality_Control/QC_Postprocessing/$case-JHU-overlay.nii.gz -l /import/monstrum/Applications/fsl/etc/luts/renderhsv.lut -z 0.5 $post/Quality_Control/QC_Postprocessing/$case-JHU-qc.png
rm -f $post/Quality_Control/QC_Postprocessing/$case-JHU-overlay.nii.gz

echo "<br> <strong>$case</strong><br>
<img src=$case-GluCEST-qc height="450" width="450">
<br>

" >> $post/Quality_Control/QC_Postprocessing/QC_GluCEST.html

echo "<br> <strong>$case</strong><br>
<img src=$case-atlas-qc.png height="400" width="400">
<br>

" >> $post/Quality_Control/QC_Postprocessing/QC_Atlases_Cortical.html

echo "<br> <strong>$case</strong><br>
<img src=$case-JHU-qc.png height="400" width="400">
<br>

" >> $post/Quality_Control/QC_Postprocessing/QC_Atlases_WM.html

#xdg-open $post/Quality_Control/QC_Postprocessing/$case-GluCEST-qc #uncomment this if you want real time GluCEST QA
#######################################################################################################

#######################################################################################################
## SVS versus CEST comparison ##

echo -e "\n Creating output files for CEST / MRS comparison for $case \n\n\n"

#Extract SVS ROI slice
extract_slice.sh -MultiLabel $SVS/$case/${case}_svsmask.nii $post/$case/orig_data/$case-B0B1CESTMAP.nii $post/$case/slices/${case}_svsmask-slice.nii

#Register SVS ROI to GluCEST space using pre-computed transform
flirt -in $post/$case/slices/${case}_svsmask-slice.nii -ref $post/$case/orig_data/$case-B0B1CESTMAP.nii -out $post/$case/$case-SVSROI-inCEST.nii.gz -init $post/$case/$case-UNI2CEST-matrix.mat -applyxfm -interp nearestneighbour -2D

#Segment SVS ROI by tissue compartment 
fslmaths $post/$case/fast/$case-2d-FAST.nii.gz -mul $post/$case/$case-SVSROI-inCEST.nii.gz $post/$case/$case-SVSROI-inCEST-seg.nii.gz

cd $post/$case
3dROIstats -mask $post/$case/$case-SVSROI-inCEST-seg.nii.gz -numROI 3 -zerofill NaN -nomeanout -nzmean -nzsigma -nzvoxels -nobriklab -1DRformat $case-CEST_b0b1thresh.nii.gz >> $post/$case/stats/$case-SVS-ROI-stats.csv

sed -i 's/name/Subject/g' $post/$case/stats/$case-SVS-ROI-stats.csv
cut -f2-3 --complement $post/$case/stats/$case-SVS-ROI-stats.csv >> tmp.csv
mv tmp.csv $post/$case/stats/$case-SVS-ROI-stats.csv

sed -i 's/\<NZMean_1\>/CSF_mean/g' $post/$case/stats/$case-SVS-ROI-stats.csv
sed -i 's/\<NZcount_1\>/CSF_numvoxels/g' $post/$case/stats/$case-SVS-ROI-stats.csv
sed -i 's/\<NZSigma_1\>/CSF_SD/g' $post/$case/stats/$case-SVS-ROI-stats.csv

sed -i 's/\<NZMean_2\>/GM_mean/g' $post/$case/stats/$case-SVS-ROI-stats.csv
sed -i 's/\<NZcount_2\>/GM_numvoxels/g' $post/$case/stats/$case-SVS-ROI-stats.csv
sed -i 's/\<NZSigma_2\>/GM_SD/g' $post/$case/stats/$case-SVS-ROI-stats.csv

sed -i 's/\<NZMean_3\>/WM_mean/g' $post/$case/stats/$case-SVS-ROI-stats.csv
sed -i 's/\<NZcount_3\>/WM_numvoxels/g' $post/$case/stats/$case-SVS-ROI-stats.csv
sed -i 's/\<NZSigma_3\>/WM_SD/g' $post/$case/stats/$case-SVS-ROI-stats.csv

if ! [ -e $post/GluCEST-SVSROI-Measures.csv ]
then
touch $post/GluCEST-SVSROI-Measures.csv
echo "Subject	CSF_mean	CSF_numvoxels	CSF_SD	GM_mean	GM_numvoxels	GM_SD	WM_mean	WM_numvoxels	WM_SD" >> $post/GluCEST-SVSROI-Measures.csv
fi

sed -n "2p" $post/$case/stats/$case-SVS-ROI-stats.csv >> $post/GluCEST-SVSROI-Measures.csv
#######################################################################################################


echo -e "\n$case SUCCESFULLY PROCESSED\n\n\n"
} | tee "$logfile"
else
echo "$case is either missing data or already processed. Will not process"
sleep 1.5
fi

done




