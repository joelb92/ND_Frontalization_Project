#!/bin/csh

#$ -q long
#$ -N frontalize_zr_morph_casia
#$ -pe smp 1
#$ -t 1-3500
set offs=0
echo $offs
set ID=`expr "$offs" + "$SGE_TASK_ID"`
echo $ID
cd /afs/crc.nd.edu/user/j/jbrogan4/Private/MATLAB/RFIIR-git/RFIIR-For_Face_Network/Vito-ViewMorphing
matlab -singleCompThread -nodisplay -nosplash -r "Frontalize_Parallel('/afs/crc.nd.edu/group/cvrl/share/jbrogan4/Datasets/PaSC_Datasets/pasc_cropped_all/test','.jpg','/scratch365/jbrogan4/pasc_webface_cmr_vitomorph','CMR','Vito',1,[250,250],'/afs/crc.nd.edu/user/j/jbrogan4/Private/MATLAB/RFIIR-git/RFIIR-For_Face_Network/Vito-ViewMorphing','pasc_cmr_morph',3500,$ID,0)"
