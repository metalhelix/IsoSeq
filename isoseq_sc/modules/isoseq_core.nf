process lima {
    tag "lima"
    input:
    val meta

    output:
    tuple path("${meta.ID}.fl.*.bam"), val(meta)

    publishDir "${meta.outFolder}/${meta.ID}/lima/", mode: 'copy'

    script:
    """
    mkdir -p ${meta.outFolder}/${meta.ID}/lima/

    ml isoseq3/4.2.0

    lima ${meta.bam} ${projectDir}/assets/primer.fasta ${meta.ID}.fl.bam --isoseq 

    # cp -f ${meta.ID}.fl.*.bam ${meta.outFolder}/${meta.ID}/lima/

    """
}

process isoSeq {
    tag "isoseq initial processing"
    errorStrategy 'retry'
    maxRetries 3
    // queue = 'bigmem'

    input:
    tuple path(fl_bam), val(meta)

    output:
    val meta 
    tuple path("${meta.ID}.dedup.*"), path("${meta.ID}.flt*") 

    publishDir "${meta.outFolder}/${meta.ID}/isoseq/", mode: 'copy'

    script:

    """
    mkdir -p ${meta.outFolder}/${meta.ID}/isoseq/ 

    ml isoseq3/4.2.0
    ml samtools

    isoseq tag ${fl_bam} ${meta.ID}.flt.bam --design T-12U-16B

    isoseq refine ${meta.ID}.flt.bam ${projectDir}/assets/primer.fasta ${meta.ID}.fltnc.bam --require-polya

    isoseq correct \
        --method percentile --percentile ${params.cells_to_filter} \
        -B ${projectDir}/assets/3M-3pgex-may-2023.REVCOMP.txt \
        ${meta.ID}.fltnc.bam \
        ${meta.ID}.fltnc.corrected.bam

    samtools sort -t CB ${meta.ID}.fltnc.corrected.bam -o ${meta.ID}.fltnc.corrected.sorted.bam

    isoseq groupdedup ${meta.ID}.fltnc.corrected.sorted.bam ${meta.ID}.dedup.bam

    # cp -f ${meta.ID}.* ${meta.outFolder}/${meta.ID}/isoseq/ 
    """
}

process pbmm2 {
    tag "pbmm2 & isoseq collapse"
    conda "/home/compbio_svc/miniconda3/envs/IsoSeq"

    input:
    val meta
    tuple path(dedup_files), path(flt_files)

    output:
    val meta
    path "${meta.ID}.*" // contain the abundance txt file
    path dedup_files

    publishDir "${meta.outFolder}/${meta.ID}/pbmm2/", mode: 'copy'

    script:
    """
    mkdir -p ${meta.outFolder}/${meta.ID}/pbmm2/

    pbmm2 align --preset ISOSEQ \
        --sort ${meta.ID}.dedup.bam \
        ${meta.ref_fasta} ${meta.ID}.ref_Mapped.bam
    
    ml isoseq3/4.2.0

    isoseq collapse ${meta.ID}.ref_Mapped.bam ${meta.ID}.collapsed.gff

    # cp ${meta.ID}.* ${meta.outFolder}/${meta.ID}/pbmm2/
    """
}

process pigeon {
    tag "pigeon"
    conda "/home/compbio_svc/miniconda3/envs/IsoSeq"

    input:
    val meta
    path all_pbmm_outputs
    path dedup_files

    output:
    val meta
    path "${meta.ID}*"
    path "seruat_dir"


    publishDir "${meta.outFolder}/${meta.ID}/pigeon/", mode: 'copy'

    script:
    """
    mkdir -p ${meta.outFolder}/${meta.ID}/pigeon/

    pigeon prepare ${meta.ID}.collapsed.gff

    pigeon classify \
        ${meta.ID}.collapsed.sorted.gff \
        ${meta.ref_gtf} \
        ${meta.ref_fasta} \
        --fl ${meta.ID}.collapsed.abundance.txt

    pigeon filter ${meta.ID}_classification.txt --isoforms ${meta.ID}.collapsed.sorted.gff --polya-percent 1 --polya-run-length 100 

    pigeon report ${meta.ID}_classification.filtered_lite_classification.txt ${meta.ID}_saturation.txt

    pigeon make-seurat --dedup ${meta.ID}.dedup.fasta --group ${meta.ID}.collapsed.group.txt -d seruat_dir ${meta.ID}_classification.filtered_lite_classification.txt

    mv seruat_dir/isoforms_seurat/genes.tsv seruat_dir/isoforms_seurat/genes.tsv.bak
    sed "s/'//g" seruat_dir/isoforms_seurat/genes.tsv.bak > seruat_dir/isoforms_seurat/genes.tsv

    mv seruat_dir/genes_seurat/genes.tsv seruat_dir/genes_seurat/genes.tsv.bak
    sed "s/'//g" seruat_dir/genes_seurat/genes.tsv.bak > seruat_dir/genes_seurat/genes.tsv

    rm seruat_dir/isoforms_seurat/genes.tsv.bak seruat_dir/genes_seurat/genes.tsv.bak
    """
} 

process qc_report {
    tag "html qc_report"
    conda "/home/compbio_svc/miniconda3/envs/SC_PRIME"

    cpus = 20
    memory { 150.GB * task.attempt }
    errorStrategy 'retry'
    maxRetries 3

    input:
    val meta
    path pigeon_allOutputs
    path seruat_dir


    output:
    path "Seurat*"

    publishDir "${meta.outFolder}/${meta.ID}/qc_report/", mode: 'copy'

    script:
    """
    mkdir -p ${meta.outFolder}/${meta.ID}/qc_report/

    cp ${projectDir}/assets/Seurat_analysis_isoform.rmd ./
    Rscript -e "rmarkdown::render('Seurat_analysis_isoform.rmd', params=list(seurat_dir='${seruat_dir}', saturation_file='${meta.ID}_saturation.txt', output_dir='${meta.outFolder}/${meta.ID}/qc_report/'))"

    cp ${projectDir}/assets/Seurat_analysis_gene.rmd ./
    Rscript -e "rmarkdown::render('Seurat_analysis_gene.rmd', params=list(seurat_dir='${seruat_dir}', saturation_file='${meta.ID}_saturation.txt', output_dir='${meta.outFolder}/${meta.ID}/qc_report/'))"
    """
}
