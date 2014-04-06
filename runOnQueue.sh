#!/bin/bash   

dataset=$1
numJob=$2
queue=$3
analysis=$4
folder=$5
where=$6
sample=$7
folder2=$8


mkdir $folder
cd $folder
cp $CMSSW_BASE/src/Analyzer/EDBRTauAnalyzer/script/createConfigFile.sh .
cp $CMSSW_BASE/src/Analyzer/EDBRTauAnalyzer/script/createRunFileCERN.sh .
cp $CMSSW_BASE/src/Analyzer/EDBRTauAnalyzer/script/createRunFileZURICH.sh .
cp $CMSSW_BASE/src/Analyzer/EDBRTauAnalyzer/script/createRunFilePERUGIA.sh .

if [ "$where" == "ZURICH" ]; then
    cp $CMSSW_BASE/src/Analyzer/EDBRTauAnalyzer/data/Zurich/$dataset\_fileList.txt .
fi
if [ "$where" == "PERUGIA" ]; then
    cp $CMSSW_BASE/src/Analyzer/EDBRTauAnalyzer/data/Perugia/$dataset\_fileList.txt .
fi

x=$(cat $dataset\_fileList.txt | wc -l)
x=$(($x+1))
maxJob=$(($x / $numJob)) 
maxJob=$(($maxJob+1))

for ((i = 0; i < $maxJob ; i++)) ;
do
    min=$(($i*$numJob));
    max=$(($(($i+1))*$numJob));
    if [ "$min" -lt "$x" ]; then
	./createConfigFile.sh $min $max $i $dataset $sample $analysis &> $dataset\_$i\_cfg.py

	if [ "$where" == "CERN" ]; then
	    ./createRunFileCERN.sh $dataset\_$i\_cfg.py $folder &> run$i.sh
	    chmod 777 run$i.sh
	    echo "bsub -q $queue run$i.sh"
	    bsub -q $queue run$i.sh
	fi
	if [ "$where" == "ZURICH" ]; then
	    ./createRunFileZURICH.sh $dataset\_$i\_cfg.py $folder &> run$i.sh
	    chmod 777 run$i.sh
	    echo "qsub -q $queue.q run$i.sh"
	    qsub -q $queue.q run$i.sh
	fi
	if [ "$where" == "PERUGIA" ]; then
	    ./createRunFilePERUGIA.sh $dataset\_$i\_cfg.py $folder $folder2 &> run$i.sh
	    chmod 777 run$i.sh
	    echo "qsub -q $queue -l nodes=sl5 run$i.sh"
	    qsub -q $queue -l nodes=sl5 run$i.sh
	fi
    fi
done
