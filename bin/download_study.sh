#!/bin/bash


study=$1

gemma-cli-sc getSingleCellDataMatrix -e $study --format mex --scale-type count --use-ensembl-ids -o $study
curl -u $GEMMA_USERNAME:$GEMMA_PASSWORD -H Accept:text/tab-separated-values https://dev.gemma.msl.ubc.ca/rest/v2/datasets/{study_name}/cellTypeAssignment?useBioAssayId=true -o $study_name.celltypes.tsv

done



