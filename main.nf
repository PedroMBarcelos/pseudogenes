nextflow.enable.dsl = 2

include { PSEUDOGENES } from './workflows/pseudogenes'

workflow {
    PSEUDOGENES()
}
