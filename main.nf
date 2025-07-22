#!/usr/bin/env nextflow

process downloadCelltypes {
    publishDir "${params.outdir}/cell_type_assignments", mode: 'copy'

    input:
       tuple val(study_name), path(study_dir)

    output:
        tuple val(study_name), path("${study_name}.celltypes.tsv"), emit: celltypes_meta

    script:
    
    def cta_protocol =  "author-submitted"
    //}
    """

   if [ ${params.author_submitted} = true ]; then

        curl -u "${params.GEMMA_USERNAME}:${params.GEMMA_PASSWORD}" \
        -H "Accept: text/tab-separated-values" \
        --compressed \
        "https://staging-gemma.msl.ubc.ca/rest/v2/datasets/${study_name}/cellTypeAssignment?useBioAssayId=true&protocol=${cta_protocol}" \
        -o "${study_name}.celltypes.tsv"
        
    else
        curl -u "${params.GEMMA_USERNAME}:${params.GEMMA_PASSWORD}" \
        -H "Accept: text/tab-separated-values" \
        --compressed \
        "https://staging-gemma.msl.ubc.ca/rest/v2/datasets/${study_name}/cellTypeAssignment?useBioAssayId=true" \
        -o "${study_name}.celltypes.tsv"
    fi
    """
}

process getGemmaMeta {
   publishDir "${params.outdir}/metadata/${study_name}", mode: 'copy'

    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
 
    input:
        tuple val(study_name), path(study_dir) 

    output:
        tuple val(study_name), path("**${study_name}_sample_meta.tsv"), emit: sample_meta

    script:

    """

    python /space/grp/rschwartz/rschwartz/get_gemma_data.nf/bin/get_gemma_meta.py \\
        --study_name ${study_name} \\

    """
}

process write_unique_cells {
    publishDir "${params.outdir}/unique_cells/${study_name}", mode: 'copy'

    input:
        tuple val(study_name), val(celltypes_meta)

    output:
        path "${study_name}_unique_cells.tsv", emit: unique_cell_path

    script:
   // metadata_file = ${meta_path}.getName()
    """
    
    awk -F'\t' '
    NR == 1 {
        for (i = 1; i <= NF; i++) if (\$i == "cell_type") col = i
        next
    }
    { count[\$col]++ }
    END {
        print "cell_type\tcount"
        for (c in count) print c "\t" count[c]
    }
    ' "$celltypes_meta" > "${study_name}_unique_cells.tsv"


    """
}

process processStudies {
    publishDir "${params.outdir}/h5ad/${study_name}", mode: 'copy'
    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
    input:
        tuple val(study_name), val(query_name), path(query_path), path(celltypes_meta), path(sample_meta)

    output:
       path "**.h5ad", emit: h5ad_paths
        
    script:

    """
    # Process the cell types metadata
    python /space/grp/rschwartz/rschwartz/get_gemma_data.nf/bin/gemma_preproc.py \\
        --query_path ${query_path} \\
        --cell_meta_path ${celltypes_meta} \\
        --sample_meta_path ${sample_meta} \\
        --query_name ${query_name} \\
        --study_name ${study_name} \\
        --gene_mapping "${params.gene_mapping}"
    """
}

include { DOWNLOAD_STUDIES_SUBWF } from "${projectDir}/modules/subworkflows/download_studies.nf"

// Workflow definition
workflow {

    
    DOWNLOAD_STUDIES_SUBWF(params.study_names, params.studies_path)

    DOWNLOAD_STUDIES_SUBWF.out.study_channel.set { study_channel }


    downloadCelltypes(study_channel)
    // Get the metadata
    getGemmaMeta(study_channel)

    
    study_channel.map { study_name, study_dir ->
        def subdirs = []
        study_dir.eachDirRecurse { it ->
            query_name = it.getName()
            query_path = it
            subdirs << [study_name, query_name, query_path]
        }
        return subdirs
    }.flatMap().set { query_paths }

   // study_dirs = downloadStudies.out.study_dir
    celltypes_meta = downloadCelltypes.out.celltypes_meta
    sample_meta = getGemmaMeta.out.sample_meta
    write_unique_cells(celltypes_meta)

    // combine all query paths with meta and cell types
    combined_all = query_paths.combine(celltypes_meta, by: 0)
    combined_all = combined_all.combine(sample_meta, by: 0)
    // Process the data
    processStudies(combined_all)
}


workflow.onError = {
println "Error: something went wrong, check the pipeline log at '.nextflow.log"
}
