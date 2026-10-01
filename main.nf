nextflow.enable.dsl = 2

include { PSEUDOGENES_TEST } from './workflows/pseudogenes_test.nf'

workflow {
    PSEUDOGENES_TEST()
}
