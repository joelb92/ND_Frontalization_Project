#!/usr/bin/env sh
# Create the imagenet lmdb inputs
# N.B. set the path to the imagenet train + val data dirs

TOOLS=/home/sandipan/caffe/.build_release/tools

TRAIN_DATA_ROOT=/home/sandipan/CW_CMR_Hass/Sym/CW_CMR_Hass_Sym_Train/
VAL_DATA_ROOT=/home/sandipan/CW_CMR_Hass/Sym/CW_CMR_Hass_Sym_Val1/
#TEST_DATA_ROOT=/home/sandipan/CW_Vito_New/Asym/CW_ZR_Val2/
# Set RESIZE=true to resize the images to 256x256. Leave as false if images have
# already been resized using another tool.
RESIZE=true
if $RESIZE; then
  RESIZE_HEIGHT=256
  RESIZE_WIDTH=256
else
  RESIZE_HEIGHT=0
  RESIZE_WIDTH=0
fi

if [ ! -d "$TRAIN_DATA_ROOT" ]; then
  echo "Error: TRAIN_DATA_ROOT is not a path to a directory: $TRAIN_DATA_ROOT"
  echo "Set the TRAIN_DATA_ROOT variable in create_imagenet.sh to the path" \
       "where the ImageNet training data is stored."
  exit 1
fi

if [ ! -d "$VAL_DATA_ROOT" ]; then
  echo "Error: VAL_DATA_ROOT is not a path to a directory: $VAL_DATA_ROOT"
  echo "Set the VAL_DATA_ROOT variable in create_imagenet.sh to the path" \
       "where the ImageNet validation data is stored."
  exit 1
fi

echo "Creating train lmdb..."

GLOG_logtostderr=1 $TOOLS/convert_imageset \
    --resize_height=$RESIZE_HEIGHT \
    --resize_width=$RESIZE_WIDTH \
    --shuffle \
	$TRAIN_DATA_ROOT \
    /home/sandipan/CW_CMR_Hass/Sym/CW_CMR_Hass_Sym_Train.txt \
    CW_cmr_hass_sym_train_lmdb

echo "Creating val1 lmdb..."

GLOG_logtostderr=1 $TOOLS/convert_imageset \
    --resize_height=$RESIZE_HEIGHT \
    --resize_width=$RESIZE_WIDTH \
    --shuffle \
	$VAL_DATA_ROOT \
    /home/sandipan/CW_CMR_Hass/Sym/CW_CMR_Hass_Sym_Val1.txt \
    CW_cmr_hass_sym_val1_lmdb
	
echo "Computing train mean..."

$TOOLS/compute_image_mean CW_cmr_hass_sym_train_lmdb CW_cmr_hass_sym_mean.binaryproto
python2.7 convert.py CW_cmr_hass_sym_mean.binaryproto CW_cmr_hass_sym_mean.npy

echo "Done."
