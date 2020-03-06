#!/bin/bash

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
    'Brassica napus' 

    
    #nice to have at some point, but missing in OMA so far:
    # 'Nicotiana tabacum' \
