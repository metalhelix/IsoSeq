include {lima; isoSeq; pbmm2; pigeon; qc_report} from './modules/isoseq_core.nf'
include {SAMPLESHEET_CHANNEL} from './workflows/samplesheet.nf' //Make sure the gff and fna are pigeon prepared! 

if (!params.samplesheet) {exit 1, "Need a samplesheet to run the workflow. Please provide a samplesheet with the --samplesheet parameter."}


workflow RUN {
    take:
    meta

    main:
    lima(meta)
    isoSeq(lima.out)
    pbmm2(isoSeq.out)
    pigeon(pbmm2.out)
    qc_report(pigeon.out)
}

workflow {
    meta = SAMPLESHEET_CHANNEL(params.samplesheet)
    RUN(meta)
}

workflow.onComplete {

    def msg = """\
        PRIME Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        workDir     : ${workflow.workDir}
        Job status : ${ workflow.success ? "Success!" : "failed, error message: ${workflow.errorMessage}" }
        """
        .stripIndent()

    println msg

    sendMail(to: 'by2747@stowers.org', subject: "Stowers ISOSEQ pipeline execution", body: msg)
}