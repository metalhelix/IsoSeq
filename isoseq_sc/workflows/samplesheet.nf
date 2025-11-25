workflow SAMPLESHEET_CHANNEL {
    take:
    samplesheet

    main:

    // Channel.from( samplesheet )
    samplesheet_ch = Channel.fromPath( samplesheet ).splitCsv( header:true, sep:',' )
        // .map{ create_metadata_channel(it) }
        // .set{ data_meta }

    emit:
    samplesheet_ch
}
