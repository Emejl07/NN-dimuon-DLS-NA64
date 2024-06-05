#!/bin/sh
# 1 parameters: 1: number of jobs, 2: jobs per events
# optional parameter: 3: config file  (to use instead of default)


#define config file to be used in this reconstruction
export CONFIGFILE=$(pwd)/config/CondorConfigSim.sh

#decide number of arguments used
if [ $# -lt 2 ] || [ $# -gt 3 ]
then
        echo "number of argument is invalid"
        echo "USAGE:"
        echo "1: NUMBER OF JOBS TO SEND TO THE CLUSTER"
        echo "2: NUMBER OF EVENTS TO SIMULATE PER JOBS"
        echo "3: CONFIG FILE TO USE INSTEAD OF DEFAULT (OPTIONAL)"
        echo "3: CURRENTLY THE CONFIG FILE USED IS:"
        echo "$CONFIGFILE"
        exit
#if number of argument is three change config file
elif [ $# -eq 3 ]
then
        CONFIGFILE=$3
        echo "MESSAGE: USING CONFIG FILE $CONFIGFILE INSTEAD OF DEFAULT"
fi

#define input
export JOBSNUMBER=$1
export EVENTSNUMBER=$2

#Take all enviromental variable from config skripts
source $CONFIGFILE

#create exepath to avoid condor incident
export EXEDIR="condor-geant4sim-$JOBNAME-$(date +%Y%m%d-%H%M%S)"
mkdir $EXEDIR

#transform condor
sed -e 's/<queue>/'$QUEUE'/g' \
        -e 's/<computing-group>/'$COMPUTINGGROUP'/g' \
        -e 's/<jobname>/'$JOBNAME'/g' \
        -e 's#<logdirectory>#'$LOGDIR'#g' \
        -e 's/<notification>/'$NOTIFICATION'/g' \
        -e 's/<usermail>/'$USERMAIL'/g' \
        -e 's/<jobsnumber>/'$JOBSNUMBER'/g' \
        $CONDORSKRIPT > "condor.submit";

#transform executable, reconstruction flag are optional
sed -e 's#<eospath>#'$EOSPATH"/"$JOBNAME'#g' \
        -e 's#<exepath>#'$LOCALPATH"/"$EXEDIR'#g' \
        -e 's#<compression>#'$COMPRESSION'#g' \
        -e 's#<compressionflags>#'$COMPRESSIONFLAGS'#g' \
        -e 's#<events>#'$EVENTSNUMBER'#g' \
        $CONDOREXE > "geant4_run.sh";

#trasform config file for geant4 simulation
sed -e 's/<verbose>/'$VERBOSE'/g' \
        -e 's/<tesla-field>/'$MAGFIELD'/g' \
        -e 's/<output-mode>/'$OUTMODE'/g' \
        -e 's/<particle>/'$PARTICLE'/g' \
        -e 's/<energy>/'$ENERGY'/g' \
        -e 's/<energy-unit>/'$ENERGYUNIT'/g' \
        -e 's/<events>/'$EVENTSNUMBER'/g' \
        $CARDEXE > "cardg.g4mac";

#give executable right
chmod a+x geant4_run.sh

#copy all needed file in the directory
cp -R -u -p $MACROPATH/mkgeant4 $LOCALPATH/$EXEDIR
cp -R -u -p $LOCALPATH/cardg.g4mac $LOCALPATH/$EXEDIR
cp -R -u -p $LOCALPATH/condor.submit $LOCALPATH/$EXEDIR
cp -R -u -p $LOCALPATH/geant4_run.sh $LOCALPATH/$EXEDIR
cp -R -u -p $SETUPFILE $LOCALPATH/$EXEDIR"/setup.sh"
cp -R -u -p $CONFIGFILE $LOCALPATH/$EXEDIR
#also copy the reconstruction in the folder if exist
#[[ -e $RECOEXEPATH/$RECOEXE ]] && cp -R -u -p $RECOEXEPATH/$RECOEXE $LOCALPATH/$EXEDIR

#call condor skript
echo "CALLING CONDOR SERVICE:"
cd $EXEDIR
condor_submit condor.submit

#back to original directory
echo "back to the original directory"
cd -
