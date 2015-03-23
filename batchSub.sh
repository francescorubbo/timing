#!/bin/bash

[ "$USER" == "pnef" ]     && WorkDir=/u/at/pnef/Work/Code/Reclustering/
[ "$USER" == "swiatlow" ] && WorkDir=/u/at/swiatlow/nfs/projects/Reclustering/
[ "$USER" == "rubbo" ]    && WorkDir=/u/at/rubbo/nfs/Timing/trunk/
# add similar line if you are not pnef

SubFileLoc=`pwd`/_batchSingleSub.sh
DateSuffix=`date +%Y%m%d_%Hh%Mmin`

echo '#!/bin/bash
echo CD to $1
echo CMD is $2

cd $1
source setup.sh
cmd=$4

echo MAKING TEMP DIR $2
JOBFILEDIR=$2
mkdir $JOBFILEDIR
REALOUT=$3
echo MADE TEMP DIR $JOBFILEDIR
echo WILL COPY TO $REALOUT

shift
shift
echo Calling $cmd $*
$cmd $*
cp -r $JOBFILEDIR/*.root $REALOUT
echo COPYING to $REALOUT
rm -rf $JOBFILEDIR
' > $SubFileLoc
chmod u+x $SubFileLoc

#----------------
Process=4
zspread=0
for mu in 80; do
    Queue=short
    nevents=200
    njobs=50
    LogPrefix=`pwd`/logs/${DateSuffix}/${DateSuffix}_bsub_${mu}_
    OutDirFinal=`pwd`/files/${DateSuffix}
    mkdir -p `dirname $LogPrefix`
    mkdir -p $OutDirFinal
    echo
    echo "Submitting $njobs jobs each with $nevents events to $Queue"
    echo $LogPrefix
    for (( ii=1; ii<=$njobs; ii++ )) ;  do
        echo $ii
        OutDir=/scratch/${DateSuffix}_${ii}/
        bsub -q ${Queue} -R rhel60 -o $LogPrefix${ii}.log $SubFileLoc           \
            ${WorkDir} ${OutDir} ${OutDirFinal} ./Timing.exe  \
            --Pileup $mu                 \
            --OutFile ${OutDir}/Sample_mu_${mu}_nevents_${nevents}_job_${ii}.root \
            --Proc ${Process} \
            --NEvents ${nevents} \
            --Zspread ${zspread} \
	    --Seed ${ii}

    done
done

