#!/bin/bash
date

dir=/star/u/gwilks3/fst/star-fwd-integration
seddir="\/star\/u\/gwilks3\/fst\/star-fwd-integration"

style=slow                 #fast, full, slow
folder=TestingMisalignment #whatever folder you used for simulation
geomtag=dev2022sm          #dev2021, dev2022, dev2022m, dev2022ms
njobs=1000                 #how many jobs you are submitting. should match njobs for submission of starsim jobs

moption=true
if [[ ${geomtag} == dev2022 || ${geomtag} == dev2021 ]]
  then
    moption=false
fi

echo "jnobs = ${njobs}"
echo "dir = ${dir}"
echo "Style = ${style}"
echo "Geometry Tag = ${geomtag}"
echo "Folder = ${folder}"

[ -d log ] || mkdir -p log
[ -d script ] || mkdir -p script
[ -d configs ] || mkdir -p configs
[ -d schedinfo ] || mkdir -p schedinfo

outdir=/gpfs01/star/pwg/gwilks3/${folder}
outtreedir=${outdir}/outtree${style}_${geomtag}
outhistdir=${outdir}/out${style}_${geomtag}
fstQAdir=${outdir}/fstQA_${style}_${geomtag}
mkdir -p ${outtreedir}
mkdir -p ${outhistdir}
mkdir -p $fstQAdir
mkdir -p $outdir/log

sed "s/geomtag/${geomtag}/g;s/dir/${seddir}/g;s/moption/${moption}/g" ${dir}/tests/${style}_track.xml > configs/${style}_track_${geomtag}_${folder}.xml         

star-submit-template -template tracking.xml -entities outdir=${outdir},outtreedir=${outtreedir},outhistdir=${outhistdir},fstQAdir=$fstQAdir,dir=$dir,geomtag=$geomtag,folder=$folder,style=$style,njobs=$njobs
