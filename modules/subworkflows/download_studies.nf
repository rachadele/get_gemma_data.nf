include { DOWNLOAD_STUDIES } from "${projectDir}/modules/processes/download_studies.nf"

workflow DOWNLOAD_STUDIES_SUBWF {

    take:
    // Exactly one of these must be non-null
    study_names   // path to a file with one study ID per line
    study_list    // comma- or space-separated study IDs passed directly on the CLI
    studies_path  // path to a pre-downloaded studies directory

    main:
    if (study_names) {
        Channel
            .fromPath(params.study_names)
            .flatMap { file -> file.readLines().collect { it.trim() }.findAll { it } }
            .set { study_names }

        DOWNLOAD_STUDIES(study_names)
            .set { study_channel }

    } else if (study_list) {
        Channel
            .from(params.study_list.split(/[,\s]+/).collect { it.trim() }.findAll { it })
            .set { study_names }

        DOWNLOAD_STUDIES(study_names)
            .set { study_channel }

    } else if (studies_path) {
        study_channel = Channel
            .fromPath(params.studies_path)
            .flatMap { path ->
                def results = []
                path.eachDir { dir -> results << [dir.name, dir] }
                return results
            }
    } else {
        exit 1, "Error: You must provide one of '--study_names', '--study_list', or '--studies_path'."
    }

    emit: study_channel
    
}
