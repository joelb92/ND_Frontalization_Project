#!/bin/csh

#$ -q long
#$ -N frontalize_zr_morph_casiaTEST
#$ -pe smp 1
#$ -t 1-100
set SGE_TASK_ID=$1
set offs=0
echo $offs
set ID=`expr "$offs" + "$SGE_TASK_ID"`
echo $ID
cd /afs/crc.nd.edu/user/j/jbrogan4/Private/MATLAB/RFIIR-git/RFIIR-For_Face_Network/Vito-ViewMorphing
matlab -nodisplay -nosplash -r "Frontalize_Parallel('/afs/crc.nd.edu/group/cvrl/share/jbrogan4/Datasets/PaSC_Datasets/pasc_cropped_all/test','.jpg','/scratch365/jbrogan4/pasc_zr_vitomorph','ZR','Vito',1,[250,250],'/afs/crc.nd.edu/user/j/jbrogan4/Private/MATLAB/RFIIR-git/RFIIR-For_Face_Network/Vito-ViewMorphing','pasc_zr_morph',5000,$ID,1)"
