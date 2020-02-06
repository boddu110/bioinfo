#!/usr/bin/env bash
#alignAll.sh
outDir='quant/'
sample='/scratch/SampleDataFiles/Paired/'
leftSuffix=".R1.paired.fastq"
rightSuffix=".R2.paired.fastq"

function align {
     for x in $sample*$leftSuffix
     do
        pathRemoved="${x/$sample/}"
        sampleName="${pathRemoved/$leftSuffix/}"
        echo $sampleName
        salmon quant -l IU \
        -1 $sample$sampleName$leftSuffix \
        -2 $sample$sampleName$rightSuffix \
        -i AipIndex \
        --validateMappings \
        -o $outDir$sampleName
        done
}

align 
#1>align.log 2>align.err &
