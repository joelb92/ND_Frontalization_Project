#!/bin/csh

#$ -q long
#$ -N frontalize_cmr_morph_casia
#$ -pe smp 1
#$ -t 1-2000
set offs=0
echo $offs
set ID=`expr "$offs" + "$SGE_TASK_ID"`
echo $ID
cd /afs/crc.nd.edu/user/j/jbrogan4/Private/MATLAB/RFIIR-git/RFIIR-For_Face_Network/Vito-ViewMorphing
matlab -nodisplay -nosplash -r "Frontalize_Parallel('/scratch365/jbrogan4/CASIA_webface_original/CASIA-WebFace','.jpg','/scratch365/jbrogan4/CASIA_webface_cmr_vitomorph','ZR','Vito',1,[250,250],'/afs/crc.nd.edu/user/j/jbrogan4/Private/MATLAB/RFIIR-git/RFIIR-For_Face_Network/Vito-ViewMorphing','CASIA_webface_cmr_morph2',2000,$ID,0)"
