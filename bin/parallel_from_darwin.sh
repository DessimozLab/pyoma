#!/bin/bash

#SBATCH --account=cdessim2_oma
#SBATCH --array=1-25
#SBATCH --partition=ax-normal
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --job-name=pyoma
#SBATCH --output=/scratch/axiom/FAC/FBM/DBC/cdessim2/oma/cur-batch/scratch/browser/logs/%A.%a.log
#SBATCH --export=None

module load oma-conv


NProc=${1:-25}
doVPS="${2:-false}"
overwrite="${3:-true}"
echo "Running with $NProc procs, doVPS=$doVPS"

if [ -z "$SLURM_ARRAY_TASK_ID" ] ; then
    echo "bla"
    SLURM_ARRAY_TASK_ID=$THIS_PROC_NR
fi
echo $SLURM_ARRAY_TASK_ID

mkdir -p ${DARWIN_NETWORK_SCRATCH_PATH}/pyoma/{prots,cps,vps}/
darwin -E -q << EOF
    NR_PROCESSES := $NProc;
    THIS_PROC_NR := $SLURM_ARRAY_TASK_ID;
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


