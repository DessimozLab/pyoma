#!/bin/bash
#SBATCH --account=cdessim2_oma
#SBATCH --partition=ax-normal
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --job-name=nebion-export
#SBATCH --output=/scratch/axiom/FAC/FBM/DBC/cdessim2/oma/cur-batch/scratch/browser/logs/nebion-%A.log
#SBATCH --export=None

module load oma-conv

python ~/pyoma/bin/export_for_genevestigator.py -v --out $DARWIN_BROWSERDATA_PATH/../downloads/nebion_orthologs.txt $DARWIN_BROWSERDATA_PATH/OmaServer.h5 \
    'Arabidopsis thaliana' \
    'Mus musculus' \
    'Rattus norvegicus' \
    'Hordeum vulgare subsp. vulgare' \
    'Homo sapiens' \
    'Macaca mulatta' \
    'Danio rerio' \
    'Oryza sativa' \
    'Glycine max' \
    'Saccharomyces cerevisiae' \
    'Zea mays' \
    'Escherichia coli' \
    'Drosophila melanogaster' \
    'Triticum aestivum' \
    'Solanum lycopersicum' \
    'Physcomitrella patens subsp. patens' \
    'Sorghum bicolor' \
    'Sus scrofa' \
    'Canis lupus familiaris' \
    'Helianthus annuus' \
    'Medicago truncatula' \
    'Brassica napus' \
    'Nicotiana tabacum'
