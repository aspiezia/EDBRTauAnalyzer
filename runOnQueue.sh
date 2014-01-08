#!/bin/bash   

dataset=$1
numJob=$2
queue=$3
analysis=$4
folder=$5
where=$6


mkdir $folder
cd $folder
cp ../script/createConfigFileFullyLep.sh .
cp ../script/createRunFileCERN.sh .
cp ../script/createRunFileZURICH.sh .
cp ../data/$dataset\_fileList.txt .

x=$(cat $dataset\_fileList.txt | wc -l)
x=$(($x+1))
maxJob=$(($x / $numJob)) 
maxJob=$(($maxJob+1))

for ((i = 0; i < $maxJob ; i++)) ;
do
    min=$(($i*$numJob));
    max=$(($(($i+1))*$numJob));
    if [ "$min" -lt "$x" ]; then
        if [ "$analysis" == "fullyLeptonic" ]; then
            ./createConfigFileFullyLep.sh $min $max $i $dataset &> $dataset\_$i\_cfg.py
        fi

	if [ "$where" = "CERN" ]; then
	    ./createRunFileCERN.sh $dataset\_$i\_cfg.py $folder &> run$i.sh
	    chmod 777 run$i.sh
	    echo "bsub -q 8nh run$i.sh"
	    bsub -q $queue run$i.sh
	fi
	if [ "$where" = "ZURICH" ]; then
	    ./createRunFileZURICH.sh $dataset\_$i\_cfg.py $folder &> run$i.sh
	    chmod 777 run$i.sh
	    echo "qsub run$i.sh"
	    qsub run$i.sh
	fi
    fi
done
