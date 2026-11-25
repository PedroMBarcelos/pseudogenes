process FILTER_TBLASTN {
    tag { "${meta.id}" }
    label 'light'

    publishDir(
        { "${params.outdir}/${meta.id}/intermediate" },
        mode: 'copy',
        pattern: '*.filtered.out'
    )

    input:
    tuple val(meta), path(tblastn_out)

    output:
    tuple val(meta),
          path("${meta.id}.tblastn.filtered.out"),
          emit: filtered_tblastn

    script:
    """
    awk -v min_identity=20 \\
        -v min_score=120 '
        NF == 0 || /^#/ {next}

        NF != 25 {
            print "Unexpected TBLASTN field count: " NF > "/dev/stderr"
            exit 2
        }

        (\$15 + 0) >= min_identity && (\$13 + 0) >= min_score {
            print
        }
    ' "${tblastn_out}" > "${meta.id}.tblastn.filtered.out"
    """
}