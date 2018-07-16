#!/bin/bash

NProc=${1:-25}
doVPS="${2:-false}"
overwrite="${3:-true}"
echo "Running with $NProc procs, doVPS=$doVPS"

if [ -z "$LSB_JOBINDEX" ] ; then
    echo "bla"
    LSB_JOBINDEX=$THIS_PROC_NR
fi
echo $LSB_JOBINDEX

mkdir -p ${DARWIN_NETWORK_SCRATCH_PATH}/pyoma/{prots,cps,vps}/
    darwin -E -q << EOF
    NR_PROCESSES := $NProc;
    THIS_PROC_NR := $LSB_JOBINDEX;
    ReadProgram('$DARWIN_OMA_REPO_PATH/lib/Platforms');
    ReadProgram('pyoma/browser/convert.drw');
    pInf := DetectParallelInfo();
    if pInf['ProcNr']=1 then
        outfn := getenv('DARWIN_NETWORK_SCRATCH_PATH').'/pyoma/gs.json';
        if $overwrite or length(FileStat(outfn))=0 then
            GetGenomeData();
        fi:
    fi:
    for g in genomes do
        if IsMyJob(pInf, g) then
            outfn := getenv('DARWIN_NETWORK_SCRATCH_PATH').'/pyoma/prots/'.g.'.json';
            if $overwrite or length(FileStat(outfn))=0 then
                GetProteinsForGenome(g); 
            fi:
            
            outfn := getenv('DARWIN_NETWORK_SCRATCH_PATH').'/pyoma/cps/'.g.'.json';
            if $overwrite or length(FileStat(outfn))=0 then
                GetSameSpeciesRelations(g); 
            fi:
            
            if $doVPS=true then 
                outfn := getenv('DARWIN_NETWORK_SCRATCH_PATH').'/pyoma/vps/'.g.'.json';
                if $overwrite or length(FileStat(outfn))=0 then
                     GetVPsForGenome(g); 
                fi:
            fi:
        fi:
    od:
    done
EOF


