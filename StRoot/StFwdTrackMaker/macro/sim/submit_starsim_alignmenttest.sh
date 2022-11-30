#!/bin/bash
date

folder=AlignmentTest_20221122_S38_noMis_noB_noV_pt5_fast_fixedKine_oldOuterCoor #whatever you want the output folder to be named
geomtag=dev2022m # dev2022, dev2022m    (ideal, misaligned)
njobs=500
nEvents=500
nTracks=1
dir=/star/u/gwilks3/fst/ForwardTracking/star-sw-1
#dir=$(echo "`pwd`" | sed 's:/:\\/:g')

echo "$dir"

logdir=/gpfs01/star/pwg/gwilks3/${folder}/log    
outdir=/gpfs01/star/pwg/gwilks3/${folder}/out 

mkdir -p ${logdir}
mkdir -p ${outdir}
mkdir -p schedinfo

star-submit-template -template ${dir}/StRoot/StFwdTrackMaker/macro/sim/starsim_alignmenttest.xml -entities nEvents=$nEvents,nTracks=$nTracks,dir=$dir,geomtag=$geomtag,njobs=$njobs,logdir=$logdir,outdir=$outdir

