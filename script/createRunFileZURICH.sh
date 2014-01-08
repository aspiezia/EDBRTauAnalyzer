#! /bin/bash

fileToRun=$1
folder=$2

echo '#!/bin/bash'
echo '#'
echo '#$ -cwd'
echo '#$ -j y'
echo '#$ -S /bin/bash'
echo '#'
echo 'date'
echo 'hostname'
echo 'pwd'
echo 'source /swshare/cms/cmsset_default.sh'
echo 'eval `scramv1 runtime -sh`'
echo "cmsRun $fileToRun"