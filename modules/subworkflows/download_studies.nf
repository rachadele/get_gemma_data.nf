include { DOWNLOAD_STUDIES } from "${projectDir}/modules/processes/download_studies.nf"

workflow DOWNLOAD_STUDIES_SUBWF {

    take:
    // Exactly one of these must be non-null
    study_names   // (--study-names) path to a file with one study ID per line
    study_file    // (--study-file)  comma- or space-separated study IDs passed directly on the CLI
    study_paths   // (--study-paths) path to a pre-downloaded studies directory

    main:
    if (study_names) {
        def study_names_ch
        if (file(study_names).exists()) {
            study_names_ch = Channel
                .fromPath(study_names)
                .flatMap { f -> f.readLines().collect { it.trim() }.findAll { it } }
        } else {
            study_names_ch = Channel
                .from(study_names.split(/[,\s]+/).collect { it.trim() }.findAll { it })
        }

        DOWNLOAD_STUDIES(study_names_ch)
            .set { study_channel }

    } else if (study_file) {
        Channel
            .from(study_file.split(/[,\s]+/).collect { it.trim() }.findAll { it })
            .set { study_names }

        DOWNLOAD_STUDIES(study_names)
            .set { study_channel }

    } else if (study_paths) {
        study_channel = Channel
            .fromPath(study_paths)
            .flatMap { path ->
                def results = []
                path.eachDir { dir -> results << [dir.name, dir] }
                return results
            }
    } else {
        exit 1, "Error: You must provide one of '--study_names', '--study_file', or '--study_paths'."
    }

    emit: study_channel
    
}
