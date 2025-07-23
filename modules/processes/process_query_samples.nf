
process PROCESS_QUERY_SAMPLE {
	tag "$query_name"
    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
    publishDir "${params.outdir}/h5ad/${study_name}", mode: 'copy', pattern: "${query_name}.h5ad"
    // add a pattern to catch small_samples/ subdirectories
    publishDir "${params.outdir}/small_samples/${study_name}", mode: 'copy', pattern: "small_samples/**"

	// This process is used to process query datasets using a pre-trained model.
	// It takes a study name and path, processes the data, and outputs processed and raw data files.
	// The model path is provided as an input parameter.


    input:
    tuple val(study_name), val(query_name), path(query_path), path(celltypes_meta), path(sample_meta)

    output:
    tuple val("${study_name}"), val("${query_name}"), path("**${query_name}.h5ad"), emit: processed_query

        
    script:


    """

    # Process the cell types metadata
    python /space/grp/rschwartz/rschwartz/get_gemma_data.nf/bin/process_query_samples.py \\
        --query_path ${query_path} \\
        --cell_meta_path ${celltypes_meta} \\
        --sample_meta_path ${sample_meta} \\
        --query_name ${query_name} \\
        --gene_mapping "${params.gene_mapping}"
    """

}