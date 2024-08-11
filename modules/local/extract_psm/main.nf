process PSMCONVERSION {
    tag "$meta.mzml_id"
    label 'process_medium'

    conda "bioconda::quantms-utils=0.0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/quantms-utils:0.0.5--pyhdfd78af_0' :
        'biocontainers/quantms-utils:0.0.5--pyhdfd78af_0' }"


    input:
    tuple val(meta), path(idxml_file), path(spectrum_df)

    output:
    path "*_psm.csv", emit: psm_info
    path "versions.yml", emit: version
    path "*.log", emit: log

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.mzml_id}"
    def string_export_decoy_psm = params.export_decoy_psm == true ? "--export_decoy_psm" : ""


    """
    quantmsutilsc psmconvert --idxml "${idxml_file}" \\
        --spectra_file ${spectrum_df} \\
        ${string_export_decoy_psm} \\
        2>&1 | tee extract_idxml.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pyopenms: \$(pip show pyopenms | grep "Version" | awk -F ': ' '{print \$2}')
    END_VERSIONS
    """
}
