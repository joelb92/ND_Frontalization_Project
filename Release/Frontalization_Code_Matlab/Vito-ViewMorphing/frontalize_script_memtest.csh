#!/bin/csh

#$ -q debug
#$ -N frontalize_zr_morph_casia_MEMTEST
#$ -pe smp 24

set SGE_TASK_ID=1003
set offs=0
echo $offs
set ID=`expr "$offs" + "$SGE_TASK_ID"`
echo $ID
cd /afs/crc.nd.edu/user/j/jbrogan4/Private/MATLAB/RFIIR-git/RFIIR-For_Face_Network/Vito-ViewMorphing
matlab -nodisplay -nosplash -r "Frontalize_Parallel('/scratch365/jbrogan4/CASIA_webface_original/CASIA-WebFace','.jpg','/scratch365/jbrogan4/CASIA_webface_zr_vitomorph','ZR','Vito',1,[250,250],'/afs/crc.nd.edu/user/j/jbrogan4/Private/MATLAB/RFIIR-git/RFIIR-For_Face_Network/Vito-ViewMorphing','CASIA_webface_zr_morph2',14000,$ID,1)"