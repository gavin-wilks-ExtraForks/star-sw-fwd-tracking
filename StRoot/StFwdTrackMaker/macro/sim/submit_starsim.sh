#!/bin/bash
date

folder=TestingFwdForAlignment_FullDetector_withB #whatever you want the output folder to be named
geomtag=dev2022m # dev2022, dev2022m    (ideal, misaligned)
njobs=2500 # each job is set for 500 events
nEvents=25
nTracks=10
dir=/star/u/gwilks3/fst/ForwardTracking/star-sw-1
#dir=$(echo "`pwd`" | sed 's:/:\\/:g')

echo "$dir"

logdir=/gpfs01/star/scratch/gwilks3/${folder}/log    
outdir=/gpfs01/star/scratch/gwilks3/${folder}/out 

mkdir -p ${logdir}
mkdir -p ${outdir}
mkdir -p schedinfo

star-submit-template -template ${dir}/StRoot/StFwdTrackMaker/macro/sim/starsim.xml -entities nEvents=$nEvents,nTracks=$nTracks,dir=$dir,geomtag=$geomtag,njobs=$njobs,logdir=$logdir,outdir=$outdir

