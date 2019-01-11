#!/bin/bash

#This script processes 7T Terra structural (UNI) data. This script assumes that relevant structural dicoms for each case exist in $dicoms, and are sorted. Furthermore, this script requires the upsampled (0.8 mm) MNI template, Harvrard Oxford atlases, JHU atlas, and Brainnetome atlas found in $path/MNI_Templates 
#This script executes the following: dcm to nifti conversion, BET masking, FAST segmentation, sub UNI to MNI registration, registration of upsampled MNI atlases to sub space using sub-specific inverse transforms

#######################################################################################################
##DEFINE PATHS##
base=/import/monstrum/SYRP/7T_Terra
path=/import/monstrum/SYRP/7T_Terra/GluCEST/Structural #path where structural output will go 
dicoms=/import/monstrum/SYRP/7T_Terra/GluCEST/Dicoms #path to dicoms
#######################################################################################################

#######################################################################################################
## IDENTIFY CASES FOR PROCESSING ##

for i in $(ls /import/monstrum/SYRP/7T_Terra/GluCEST/Dicoms/) # script will only be executed for cases in this directory
do
case=${i##*/}
echo "CASE: $case"

if ! [ -d $path/$case ] && [ -d $dicoms/$case/*mp2rage_mark_ipat3_0.80mm_INV2 ] && [ -d $dicoms/$case/*mp2rage_mark_ipat3_0.80mm_UNI_Images ]

then
echo "--------Processing structural data for $case---------" 
sleep 1.5
mkdir $path/$case
mkdir $path/$case/fast
mkdir $path/$case/MNI_transforms
mkdir $path/$case/atlases
#######################################################################################################

#######################################################################################################
echo "----------Converting dicoms for $case----------"
#convert UNI
dcm2nii -d N -e N -f N -i N -r N -x N -o $path/$case $dicoms/$case/*mp2rage_mark_ipat3_0.80mm_UNI_Images
mv $path/$case/mp2ragemarkipat3080mm.nii.gz $path/$case/$case-UNI.nii.gz

#convert INV2 for mask generatation
dcm2nii -d N -e N -f N -i N -r N -x N -o $path/$case $dicoms/$case/*mp2rage_mark_ipat3_0.80mm_INV2
mv $path/$case/mp2ragemarkipat3080mm.nii.gz $path/$case/$case-INV2.nii.gz
#######################################################################################################

#######################################################################################################
echo "----------Masking $case----------"
#mask INV2 with BET
bet $path/$case/$case-INV2.nii.gz $path/$case/$case -m -f 0.2 #using INV2 as input and -f 0.2 determined by David and Eric after testing BET masking on INV2, INV1 and UNI images and testing a range of options for -f
rm -f $path/$case/$case.nii.gz
fslmaths $path/$case/$case-UNI.nii.gz -mul $path/$case/${case}_mask.nii.gz $path/$case/$case-UNI.nii.gz
rm -f $path/$case/${case}_mask.nii.gz $path/$case/$case-INV2.nii.gz

#Create and erode brain mask
fslmaths $path/$case/$case-UNI.nii.gz -bin $path/$case/$case-mask.nii.gz
fslmaths $path/$case/$case-mask.nii.gz -ero -kernel sphere 1 $path/$case/$case-UNI-mask-er.nii.gz
rm -f $path/$case/$case-mask.nii.gz 

#Apply new eroded mask to UNI image
fslmaths $path/$case/$case-UNI.nii.gz -mas $path/$case/$case-UNI-mask-er.nii.gz $path/$case/$case-UNI-masked.nii.gz
#######################################################################################################

#######################################################################################################
echo "----------Performing FAST on $case----------"
fast -n 3 -t 1 -g -o $path/$case/fast/$case $path/$case/$case-UNI-masked.nii.gz
#######################################################################################################

#######################################################################################################
echo "----------Registering $case to MNI space----------"
#Register masked UNI to upsampled MNI T1 template
#Note: MNI 152 T1 template was upsampled as follows: ResampleImage 3 MNI152lin_T1_1mm_brain.nii.gz MNI152lin_T1_0.8mm_brain.nii.gz 0.8223684430X0.8223684430X0.8199999928 0 4

antsRegistrationSyN.sh -d 3 -f $path/MNI_Templates/MNI/MNI152lin_T1_0.8mm_brain.nii.gz -m $path/$case/$case-UNI-masked.nii.gz -o $path/$case/MNI_transforms/$case-UNIinMNI-

echo "----------Applying Transforms to $case to generate atlases----------"
#Apply transforms to Brainnetome atlas
#Note: upsampled atlas used: ResampleImage 3 BNA-maxprob-thr25-1mm.nii.gz BNA-maxprob-thr25-0.8mm 0.8223684430X0.8223684430X0.8199999928 0 1

antsApplyTransforms -d 3 -r $path/$case/$case-UNI-masked.nii.gz -i $path/MNI_Templates/Brainnetome/BNA-maxprob-thr25-0.8mm.nii.gz -n MultiLabel -o $path/$case/atlases/$case-BNA-tmp.nii.gz -t [$path/$case/MNI_transforms/$case-UNIinMNI-0GenericAffine.mat,1] -t $path/$case/MNI_transforms/$case-UNIinMNI-1InverseWarp.nii.gz

#Apply transforms to HarvardOxford atlases
#Note: upsampled atlas used, upsampled as above ^^

	for reg in cort sub
	do
	antsApplyTransforms -d 3 -r $path/$case/$case-UNI-masked.nii.gz -i $path/MNI_Templates/HarvardOxford/HarvardOxford-$reg-maxprob-thr25-0.8mm.nii.gz -n MultiLabel -o $path/$case/atlases/$case-HarvardOxford-$reg-tmp.nii.gz -t [$path/$case/MNI_transforms/$case-UNIinMNI-0GenericAffine.mat,1] -t $path/$case/MNI_transforms/$case-UNIinMNI-1InverseWarp.nii.gz
	done

echo "----------Multiplying FAST GM segmentation and atlas for $case----------"
#Increase the precision of the gray matter parcellations 

fslmaths $path/$case/atlases/$case-BNA-tmp.nii.gz -mul $path/$case/fast/${case}_seg_1.nii.gz  $path/$case/atlases/$case-BNA.nii.gz

	for reg in cort sub
	do
	fslmaths $path/$case/atlases/$case-HarvardOxford-$reg-tmp.nii.gz -mul $path/$case/fast/${case}_seg_1.nii.gz $path/$case/atlases/$case-HarvardOxford-$reg.nii.gz
	done

echo "----------Applying Transforms to $case to generate WM atlas----------"
#atlas upsampled with ResampleImage 3 JHU-ICBM-labels-1mm.nii.gz JHU-ICBM-labels-0.8mm.nii.gz 0.8223684430X0.8223684430X0.8199999928 0 1
antsApplyTransforms -d 3 -r $path/$case/$case-UNI-masked.nii.gz -i $path/MNI_Templates/JHU/JHU-ICBM-labels-0.8mm.nii.gz -n MultiLabel -o $path/$case/atlases/$case-JHU-tmp.nii.gz -t [$path/$case/MNI_transforms/$case-UNIinMNI-0GenericAffine.mat,1] -t $path/$case/MNI_transforms/$case-UNIinMNI-1InverseWarp.nii.gz


echo "----------Multiplying FAST WM segmentation and atlas for $case----------"
#Increase the precision of the white matter parcellations 

fslmaths $path/$case/atlases/$case-JHU-tmp.nii.gz -mul $path/$case/fast/${case}_seg_2.nii.gz  $path/$case/atlases/$case-JHU.nii.gz
#######################################################################################################

#clean up
rm -f $path/$case/atlases/$case-JHU-tmp.nii.gz


#clean up
rm -f $path/$case/$case-UNI.nii.gz
rm -f $path/$case/atlases/*tmp.nii.gz

echo -e "\n$case SUCCESFULLY PROCESSED\n\n\n"

else
echo "$case is either missing structural dicoms or already processed. Will not process"
sleep 1.5
fi

done 




