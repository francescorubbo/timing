#!/bin/bash

[ "$USER" == "kurinsky" ] && WorkDir=/u/ki/kurinsky/ATLAS/
# add similar line if you are not kurinsky

SubFileLoc=`pwd`/_submitSingleJob.sh
DateSuffix=`date +%Y%m%d_%Hh%Mmin`

export BSUB_QUIET=

echo '#!/bin/bash
echo CD to $1
echo CMD is $4

cd $1
cmd=$4

echo MAKING TEMP DIR $2
JOBFILEDIR=$2
mkdir $JOBFILEDIR
REALOUT=$3
echo MADE TEMP DIR $JOBFILEDIR
echo WILL COPY TO $REALOUT

shift 4
echo Calling $cmd $*
$cmd $*
cp -r $JOBFILEDIR/*.root $REALOUT
echo COPYING to $REALOUT
rm -rf $JOBFILEDIR
' > $SubFileLoc
chmod u+x $SubFileLoc

#----------------
# OPTIONS

Process=4
bunchsize="0.075"
psis="5"
mus="200"
pixelSizes="500"
profiles="1"
timeModes="0"
Queue=short
nevents=20
njobs=500
HSMode="SmearHSZT"
PUMode="VaryZT"
flags="--TrueVelocity"

OutDirFinal=`pwd`/files/${DateSuffix}
mkdir -p $OutDirFinal
echo
echo "Submitting $njobs jobs each with $nevents events to $Queue"

echo "Job Parameter Grid:"
echo "Psi - "$psis
echo "Mu  - "$mus
echo "PixelSize - "$pixelSizes
echo "Profile - "$profiles
echo "JetTiming - "$timeModes
echo "Bunchsize - "$bunchsize
echo "HSMode - "$HSMode
echo "PUMode - "$PUMode
echo "Extra Flags - "$flags

for psi in $psis ; do
    for mu in $mus ; do
	for profile in $profiles ; do
	    for pixelSize in $pixelSizes ; do
		for tm in $timeModes ; do
		    LogPrefix=`pwd`/logs/${DateSuffix}/${DateSuffix}_${mu}_${psi}_${nevents}_${profile}_px-${pixelSize}_time-${tm}
		    mkdir -p `dirname $LogPrefix`
		    echo $LogPrefix
		    
		    for (( ii=1; ii<=$njobs; ii++ )) ;  do
			OutDir=/scratch/${DateSuffix}_${ii}/
			
			bsub -q ${Queue} -R rhel60 -o $LogPrefix${ii}.log \
			    $SubFileLoc ${WorkDir} ${OutDir} ${OutDirFinal} \
			    Timing.sh  \
			    --Pileup $mu                 \
			    --OutFile ${OutDir}/Sample_mu-${mu}_psi-${psi}_nevents-${nevents}_prof-${profile}_px-${pixelSize}_time-${tm}_job-${ii}.root \
			    --Proc ${Process} \
			    --NEvents ${nevents} \
			    --BunchSize ${bunchsize} \
			    --Seed ${ii} \
			    --${HSMode} \
			    --${PUMode} \
			    --Profile $profile \
			    --Psi $psi   \
			    --PixelSize $pixelSize \
			    --JetTiming ${tm} \
			    $flags
		    done
		done
	    done
	done
    done
done