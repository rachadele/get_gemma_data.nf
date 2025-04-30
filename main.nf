#!/usr/bin/env nextflow




 process downloadStudies {
    //publishDir "${params.outdir}/studies", mode: 'copy'

    input:
        val study_name

    output:
        tuple val(study_name), path("${study_name}/"), emit: study_dir


     script:

     """
     gemma-cli-sc getSingleCellDataMatrix --no-streaming -e $study_name --format mex --scale-type count --use-ensembl-ids -o $study_name 


     """
 }

process downloadCelltypes {
    publishDir "${params.outdir}/cell_type_assignments", mode: 'copy'

    input:
        val study_name

    output:
        tuple val(study_name), path("${study_name}.celltypes.tsv"), emit: celltypes_meta

    script:
    
    """
    # if params.author_submitted is true, use the author-submitted cell type assignments
    # if params.author_submitted is false, use the preferred cell type assignments from the single cell dimension

   if [ ${params.author_submitted} = true ]; then

        curl -u "${params.GEMMA_USERNAME}:${params.GEMMA_PASSWORD}" \
        -H "Accept: text/tab-separated-values" \
        --compressed \
        "https://dev.gemma.msl.ubc.ca/rest/v2/datasets/${study_name}/cellTypeAssignment?useBioAssayId=true&cellTypeAssignment=author-submitted" \
        -o "${study_name}.celltypes.tsv"

    else
        curl -u "${params.GEMMA_USERNAME}:${params.GEMMA_PASSWORD}" \
        -H "Accept: text/tab-separated-values" \
        --compressed \
        "https://dev.gemma.msl.ubc.ca/rest/v2/datasets/${study_name}/cellTypeAssignment?useBioAssayId=true" \
        -o "${study_name}.celltypes.tsv"
    fi
    """

    //        #curl -u "$GEMMA_USERNAME:$GEMMA_PASSWORD" -H "Accept:text/tab-separated-values" \
    // "https://dev.gemma.msl.ubc.ca/rest/v2/datasets/GSE198014/singleCellDimension" \
   //  -o GSE198014.celltypes.tsv
   // awk -F'\t' '$2 == 2' GSE198014.celltypes.tsv > GSE198014.celltypes_filtered.tsv
}

process getGemmaMeta {
    publishDir "${params.outdir}", mode: 'copy'

    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
 
    input:
        val study_name

    output:
        tuple val(study_name), path("**${study_name}_sample_meta.tsv"), emit: sample_meta

    script:

    """

    python /space/grp/rschwartz/rschwartz/get_gemma_data.nf/bin/get_gemma_meta.py \\
        --study_name ${study_name} \\

    """
}


process processStudies {
    publishDir "${params.outdir}", mode: 'copy'
    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
    input:
        tuple val(study_name), path(study_dir), path(celltypes_meta), path(sample_meta)

    output:
       path "**.h5ad", emit: h5ad_paths
       path "**${study_name}_unique_cells.tsv", emit: unique_cell_path
        
    script:

    """
    # Process the cell types metadata
    python /space/grp/rschwartz/rschwartz/get_gemma_data.nf/bin/gemma_preprop.py \\
        --study_dir ${study_dir} \\
        --cell_meta_path ${celltypes_meta} \\
        --sample_meta_path ${sample_meta} \\
        --study_name ${study_name} \\
        ${params.write_samples ? "--write_samples" : ""}
    """
}




// Workflow definition
workflow {

    // Define the study names
    study_names = Channel.fromPath(params.study_names).flatMap { file ->
        // Read the file, split by lines, and trim any extra spaces
        file.readLines().collect { it.trim() }
    }


    // Download the data
    downloadStudies(study_names)
    downloadCelltypes(study_names)
    // Get the metadata
    getGemmaMeta(study_names)

    study_dirs = downloadStudies.out.study_dir
    celltypes_meta = downloadCelltypes.out.celltypes_meta
    sample_meta = getGemmaMeta.out.sample_meta

    combined_params_1 = study_dirs.combine(celltypes_meta, by: 0)
    combined_params_2 = combined_params_1.combine(sample_meta, by: 0)
    //combine all of these based on study name


    // Process the data
    processStudies(combined_params_2)

}


workflow.onError = {
println "Error: something went wrong, check the pipeline log at '.nextflow.log"
}
