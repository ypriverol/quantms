// OpenMS Unimod name -> monoisotopic delta mass (Da)
def ANDES_MOD_MASS = [
    'Carbamidomethyl': 57.02146,
    'Oxidation'      : 15.99491,
    'Acetyl'         : 42.01057,
    'Phospho'        : 79.96633,
    'TMT6plex'       : 229.16293,
    'TMT10plex'      : 229.16293,
    'TMT16plex'      : 304.20715,
    'TMT18plex'      : 304.20715,
    'iTRAQ4plex'     : 144.10206,
    'iTRAQ8plex'     : 304.20536,
    'Deamidated'     : 0.98402,
]

// "Oxidation (M)" / "Acetyl (Protein N-term)" / "TMT6plex (N-term)" -> andes mods.txt line
def andesModLine = { String openmsMod, String kind ->   // kind = 'fix' | 'opt'
    def m = (openmsMod =~ /^(.+?)\s*\((.+)\)\s*$/)
    assert m.matches() : "Unrecognised modification: ${openmsMod}"
    def name = m.group(1).trim()
    def target = m.group(2).trim()
    def mass = ANDES_MOD_MASS[name]
    assert mass != null : "No mass mapping for modification '${name}'. Add it to ANDES_MOD_MASS."
    def residue = '*'
    def location = 'any'
    switch (target.toLowerCase()) {
        case 'protein n-term': residue='*'; location='Prot-N-term'; break
        case 'protein c-term': residue='*'; location='Prot-C-term'; break
        case 'n-term':         residue='*'; location='N-term';      break
        case 'c-term':         residue='*'; location='C-term';      break
        default:               residue=target; location='any';      break  // single residue e.g. "M","C","K"
    }
    return "${mass},${residue},${kind},${location},${name}"
}

process ANDES {
    tag "$meta.mzml_id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        ( params.andes_container ?: 'oras://ghcr.io/bigbio/andes-sif:0.1.0' ) :
        ( params.andes_container ?: 'ghcr.io/bigbio/andes:0.1.0' ) }"

    input:
    tuple val(meta), path(mzml_file), path(database)

    output:
    tuple val(meta), path("${mzml_file.baseName}_andes.idparquet"), emit: id_files_andes
    path "versions.yml", emit: versions
    path "*.log",        emit: log
    path "mods.txt",     emit: mods_file

    script:
    def args = task.ext.args ?: ''

    // enzyme name -> andes enzyme slug
    def enzymeMap = [ 'Trypsin':'trypsin', 'Trypsin/P':'trypsin', 'Arg-C':'argc',
                      'Asp-N':'aspn', 'Chymotrypsin':'chymotrypsin', 'Lys-C':'lysc',
                      'Lys-N':'lysn', 'Glu-C':'gluc', 'unspecific cleavage':'nonspecific' ]
    def andesEnzyme = enzymeMap[meta.enzyme] ?: 'trypsin'

    def ntt = (meta.enzyme == 'unspecific cleavage') ? 'non-specific' :
              (params.num_enzyme_termini == 'fully') ? 'fully' :
              (params.num_enzyme_termini == 'none')  ? 'non-specific' : 'semi'

    def iso = params.isotope_error_range.split(',')
    def isoMin = iso[0].trim(); def isoMax = iso[1].trim()

    // precursor tolerance: ppm vs Da
    def precFlag = (meta.precursormasstoleranceunit == 'ppm') ?
        "--precursor-tol-ppm ${meta.precursormasstolerance}" :
        "--precursor-tol-da ${meta.precursormasstolerance}"

    // build mods.txt from meta (fixed + variable)
    def fixedLines = meta.fixedmodifications?.trim() ?
        meta.fixedmodifications.tokenize(',').collect { andesModLine(it.trim(), 'fix') } : []
    def varLines = meta.variablemodifications?.trim() ?
        meta.variablemodifications.tokenize(',').collect { andesModLine(it.trim(), 'opt') } : []
    def modsContent = (["NumMods=${params.max_mods}"] + fixedLines + varLines).join('\n')

    def scoreFlag    = params.andes_score == 'strong' ? '--score strong' : '--score rank'
    def chimericFlag = params.andes_chimeric ? '--chimeric' : ''
    def refineFlag   = params.andes_refine   ? '--refine'   : ''
    def optFlags     = [scoreFlag, chimericFlag, refineFlag].findAll { it }.join(' ')

    """
    cat > mods.txt <<'EOF'
${modsContent}
EOF

    andes \\
        --spectrum ${mzml_file} \\
        --database "${database}" \\
        --output-pin ${mzml_file.baseName}_andes.pin \\
        --output-parquet ${mzml_file.baseName}_andes.idparquet \\
        --threads $task.cpus \\
        --enzyme ${andesEnzyme} \\
        --enzyme-specificity ${ntt} \\
        --max-missed-cleavages $params.allowed_missed_cleavages \\
        --min-length $params.min_peptide_length \\
        --max-length $params.max_peptide_length \\
        --top-n $params.num_hits \\
        --max-mods $params.max_mods \\
        --charge-min $params.min_precursor_charge \\
        --charge-max $params.max_precursor_charge \\
        --isotope-error-min ${isoMin} \\
        --isotope-error-max ${isoMax} \\
        ${precFlag} \\
        --mods mods.txt \\
        ${optFlags} \\
        $args \\
        2>&1 | tee ${mzml_file.baseName}_andes.log

    # andes has no --version flag yet; pin until the binary exposes one
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        andes: 0.1.0
    END_VERSIONS
    """
}
