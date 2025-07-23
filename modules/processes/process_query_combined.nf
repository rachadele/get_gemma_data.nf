
process PROCESS_QUERY_COMBINED {
    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
	tag "$study_name"

	// This process is used to process query datasets using a pre-trained model.
	// It takes a study name and path, processes the data, and outputs processed and raw data files.
	// The model path is provided as an input parameter.


    input:
    tuple val(study_name), path(study_path), path(celltypes_meta), path(sample_meta)

    output:
    tuple val("${study_name}"), val("${study_name}"), path("${study_name}.h5ad"), emit: processed_query

        
    script:


    """

    python $projectDir/bin/process_query.py \\
                --study_path ${study_path} \\
                --study_name ${study_name} \\
                --cell_meta_path ${celltypes_meta} \\
                --sample_meta_path ${sample_meta} \\
                --gene_mapping "${params.gene_mapping}" 

     """

}