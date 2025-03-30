//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA                                                 } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_NCRNA_FASTA                                           } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF                                                   } from '../../modules/nf-core/gunzip/main'
include { UNTAR as UNTAR_BT                                                      } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_STAR                                                    } from '../../modules/nf-core/untar/main'
include { BOWTIE_BUILD                                                           } from '../../modules/nf-core/bowtie/build/main'
include { STAR_GENOMEGENERATE                                                    } from '../../modules/nf-core/star/genomegenerate/main'
include { SAMTOOLS_FAIDX as GENOME_INDEX                                         } from '../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FAIDX as NCRNA_INDEX                                          } from '../../modules/nf-core/samtools/faidx/main'
include { LINUX_COMMAND as REMOVE_GTF_BRACKETS                                   } from '../../modules/local/linux_command'
include { CUSTOM_GETCHROMSIZES as GENOME_CHROM_SIZE                              } from '../../modules/nf-core/custom/getchromsizes/main'
include { CUSTOM_GETCHROMSIZES as NCRNA_CHROM_SIZE                               } from '../../modules/nf-core/custom/getchromsizes/main'
include { PREPARE_REPRESENTATIVE_TRANSCRIPT                                      } from '../../modules/local/find_longest_transcript/main'
include { FILTER_GTF_BY_TRANSCRIPT                                               } from '../../modules/local/filter_gtf/main'
include { ICOUNTMINI_SEGMENT as ICOUNT_SEG_GTF                                   } from '../../modules/nf-core/icountmini/segment/main'
include { ICOUNTMINI_SEGMENT as ICOUNT_SEG_FILTGTF                               } from '../../modules/nf-core/icountmini/segment/main'
include { CLIPSEQ_RESOLVE_UNANNOTATED as RESOLVE_UNANNOTATED_REGIONS             } from '../../modules/local/resolve_unannotated/main'

workflow PREPARE_GENOME {
    take:
    fasta                      // file: .fasta
    fasta_fai                  // file: .fai
    ncrna_fasta                // file: .fasta
    ncrna_fasta_fai            // file: .fai
    gtf                        // file: .gtf
    genome_index               // folder: index
    ncrna_genome_index         // folder: index
    genome_chrom_sizes         // file: .txt
    ncrna_chrom_sizes          // file: .txt
    representative_transcript  // file: .txt
    representative_transcript_fai  // file: .fai
    representative_transcript_gtf  // file: .gtf
    filtered_gtf               // file: .gtf
    seg_gtf                    // file: .gtf
    regions_gtf                // file: .gtf
    regions_filt_gtf           // file: .gtf
    regions_resolved_gtf       // file: .gtf

    main:

    // Init
    ch_versions = Channel.empty()
    //
    // MODULE: Uncompress genome fasta file if required
    //
    ch_fasta = Channel.empty()
    if (fasta.toString().endsWith(".gz")) {
        ch_fasta    = GUNZIP_FASTA ( [ [id:fasta.baseName], fasta ] ).gunzip
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.of([ [id:fasta.baseName], fasta ])
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], fasta]
    //ch_fasta | view

    //
    // MODULE: Uncompress genome ncrna_fasta file if required
    //
    ch_ncrna_fasta = Channel.empty()
    if (ncrna_fasta.toString().endsWith(".gz")) {
        ch_ncrna_fasta = GUNZIP_NCRNA_FASTA ( [ [id:ncrna_fasta.baseName], ncrna_fasta ] ).gunzip
        ch_versions = ch_versions.mix(GUNZIP_NCRNA_FASTA.out.versions)
    } else {
        ch_ncrna_fasta = Channel.of([ [id:ncrna_fasta.baseName], ncrna_fasta ])
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], fasta]
    //ch_ncrna_fasta | view

    //
    // MODULE: Uncompress genome gtf file if required
    //
    ch_gtf = Channel.empty()
    if (gtf.toString().endsWith(".gz")) {
        ch_gtf      = GUNZIP_GTF ( [ [id:gtf.baseName], gtf ] ).gunzip
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    } else {
        ch_gtf = Channel.of([ [id:gtf.baseName], gtf ])
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], gtf]
    //ch_gtf | view

    //
    // MODULES: Uncompress STAR index or generate if required
    //
    ch_star_index = Channel.empty()
    if (genome_index) {
        if (genome_index.toString().endsWith(".tar.gz")) {
            ch_star_index = UNTAR_STAR ( [ [:], genome_index ] ).untar
            ch_versions  = ch_versions.mix(UNTAR_STAR.out.versions)
        } else {
            ch_star_index = Channel.of([ [:] , genome_index ])
        }
    }
    else {
        ch_star_index = STAR_GENOMEGENERATE ( ch_fasta, ch_gtf ).index
        ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
    }

    //
    // MODULES: Uncompress Bowtie index or generate if required
    //
    ch_bt_index = Channel.empty()
    if (ncrna_genome_index) {
        if (ncrna_genome_index.toString().endsWith(".tar.gz")) {
            ch_bt_index = UNTAR_BT ( [ [:], ncrna_genome_index ] ).untar
            ch_versions  = ch_versions.mix(UNTAR_BT.out.versions)
        } else {
            ch_bt_index = Channel.of([ [:] , ncrna_genome_index ])
        }
    }
    else {
        ch_bt_index = BOWTIE_BUILD ( ch_ncrna_fasta ).index
        ch_versions = ch_versions.mix(BOWTIE_BUILD.out.versions)
    }

    //
    // MODULE: Create fasta fai if required
    //
    ch_fasta_fai = Channel.empty()
    if (fasta_fai) {
        ch_fasta_fai = Channel.of([ [id:fasta_fai.baseName], fasta_fai ])
    } else {
        GENOME_INDEX (
            ch_fasta,
            [[],[]]
        )
        ch_fasta_fai = GENOME_INDEX.out.fai
        ch_versions = ch_versions.mix(GENOME_INDEX.out.versions)
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], fai]
    //ch_fasta_fai | view

    //
    // MODULE: Create fasta fai if required for ncrna genome
    //
    ch_ncrna_fasta_fai = Channel.empty()
    if (fasta_fai) {
        ch_ncrna_fasta_fai = Channel.of([ [id:ncrna_fasta_fai.baseName], fasta_fai ])
    } else {
        NCRNA_INDEX (
            ch_ncrna_fasta,
            [[],[]]
        )
        ch_ncrna_fasta_fai = NCRNA_INDEX.out.fai
        ch_versions = ch_versions.mix(NCRNA_INDEX.out.versions)
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], fai]
    //ch_fasta_fai | view

    //
    // MODULE: Calc genome chrom sizes
    //
    ch_genome_chrom_sizes = genome_chrom_sizes
    if(!genome_chrom_sizes) {
        ch_genome_chrom_sizes = GENOME_CHROM_SIZE ( ch_fasta ).sizes
        ch_versions  = ch_versions.mix(GENOME_CHROM_SIZE.out.versions)
    }

    //
    // MODULE: Calc ncrna chrom sizes
    //
    ch_ncrna_chrom_sizes = ncrna_chrom_sizes
    if(!ncrna_chrom_sizes) {
        ch_ncrna_chrom_sizes = NCRNA_CHROM_SIZE ( ch_ncrna_fasta ).sizes
        ch_versions  = ch_versions.mix(NCRNA_CHROM_SIZE.out.versions)
    }

    //
    // MODULE: Remove brackets from in gene names from GTF as causes UMICollapse to fail.
    //
    REMOVE_GTF_BRACKETS (
        ch_gtf,
        [],
        false
    )
    ch_gtf = REMOVE_GTF_BRACKETS.out.file
    // EXAMPLE CHANNEL STRUCT: [[meta], fai]
    //ch_gtf | view

    //
    // MODULE: Validate provided transcripts and generate transcript files (fai, GTF) if needed. Uses provided transcripts or selects longest from primary genome if missing.
    //

    // Channel for representative_transcript
    ch_representative_transcript = representative_transcript ?
        Channel.of([ [id: representative_transcript.baseName], representative_transcript ]) :
        Channel.of([[], []])

    // Similarly for fai and gtf, but set as empty if not provided
    ch_representative_transcript_fai = representative_transcript_fai ?
        Channel.of([ [id: representative_transcript_fai.baseName], representative_transcript_fai ]) :
        Channel.empty()

    ch_representative_transcript_gtf = representative_transcript_gtf ?
        Channel.of([ [id: representative_transcript_gtf.baseName], representative_transcript_gtf ]) :
        Channel.empty()

    // Run PREPARE_REPRESENTATIVE_TRANSCRIPT: this process is always executed (even if files given) to validate provided representative transcripts and ensure consistency with the genome GTF
    PREPARE_REPRESENTATIVE_TRANSCRIPT (
        ch_gtf,
        ch_representative_transcript
    )

    ch_representative_transcript     = PREPARE_REPRESENTATIVE_TRANSCRIPT.out.representative_transcript
    ch_representative_transcript_fai = PREPARE_REPRESENTATIVE_TRANSCRIPT.out.representative_transcript_fai
    ch_representative_transcript_gtf = PREPARE_REPRESENTATIVE_TRANSCRIPT.out.representative_transcript_gtf
    ch_versions                      = ch_versions.mix(PREPARE_REPRESENTATIVE_TRANSCRIPT.out.versions)

    //
    // MODULE: Filter GTF by representative transcripts
    //
    ch_filt_gtf = Channel.of( [ [id:filtered_gtf.baseName], filtered_gtf ] )
    if (!filtered_gtf) {
        FILTER_GTF_BY_TRANSCRIPT (
            ch_gtf,
            ch_representative_transcript
        )
        ch_filt_gtf = FILTER_GTF_BY_TRANSCRIPT.out.filtered_gtf
        ch_versions = ch_versions.mix(FILTER_GTF_BY_TRANSCRIPT.out.versions)
    }

    //
    // MODULE: Segment GTF file using icount
    //
    ch_seg_gtf     = Channel.of( [ [id:seg_gtf.baseName], seg_gtf ] )
    ch_regions_gtf = Channel.of( [ [id:regions_gtf.baseName], regions_gtf ] )
    if (!seg_gtf || !regions_gtf) {
        ICOUNT_SEG_GTF (
            ch_gtf,
            ch_fasta_fai.map{ it[1] }
        )
        ch_seg_gtf     = ICOUNT_SEG_GTF.out.gtf
        ch_regions_gtf = ICOUNT_SEG_GTF.out.regions
        ch_versions    = ch_versions.mix(ICOUNT_SEG_GTF.out.versions)
    }

    // EXAMPLE CHANNEL STRUCT: [[meta], gtf]
    //ICOUNT_SEG_GTF.out.gtf | view

    //
    // MODULE: Segment the filtered GTF file using icount
    //
    ch_regions_filt_gtf = Channel.of( [ [id:regions_filt_gtf.baseName], regions_filt_gtf ] )
    if (!regions_filt_gtf) {
        ICOUNT_SEG_FILTGTF (
            ch_filt_gtf,
            ch_fasta_fai.map{ it[1] }
        )
        ch_regions_filt_gtf = ICOUNT_SEG_FILTGTF.out.regions
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], gtf]
    //ICOUNT_SEG_FILTGTF.out.gtf | view


    //
    // MODULE: Resolve the GTF regions that iCount did not annotate REGIONS FILE
    //
    ch_regions_resolved_gtf = Channel.of( [ [id:regions_resolved_gtf.baseName], regions_resolved_gtf ] )
    if (!regions_resolved_gtf) {
        RESOLVE_UNANNOTATED_REGIONS (
            ch_regions_gtf,
            ch_regions_filt_gtf,
            ch_fasta_fai
        )
        ch_regions_resolved_gtf = RESOLVE_UNANNOTATED_REGIONS.out.gtf
    }
    // EXAMPLE CHANNEL STRUCT: [[meta], gtf]
    //RESOLVE_UNANNOTATED_REGIONS.out.gtf | view



    emit:
    fasta                         = ch_fasta                             // channel: [ val(meta), [ fasta ] ]
    fasta_fai                     = ch_fasta_fai                         // channel: [ val(meta), [ fai ] ]
    ncrna_fasta                   = ch_ncrna_fasta                       // channel: [ val(meta), [ fasta ] ]
    ncrna_fasta_fai               = ch_ncrna_fasta_fai                   // channel: [ val(meta), [ fai ] ]
    genome_index                  = ch_star_index                        // channel: [ val(meta), [ star_index ] ]
    ncrna_index                   = ch_bt_index                          // channel: [ val(meta), [ bt2_index ] ]
    chrom_sizes                   = ch_genome_chrom_sizes                // channel: [ val(meta), [ txt ] ]
    ncrna_chrom_sizes             = ch_ncrna_chrom_sizes                 // channel: [ val(meta), [ txt ] ]
    gtf                           = ch_gtf                               // channel: [ val(meta), [ gtf ] ]
    representative_transcript     = ch_representative_transcript         // channel: [ val(meta), [ txt ] ]
    representative_transcript_fai = ch_representative_transcript_fai     // channel: [ val(meta), [ fai ] ]
    representative_transcript_gtf = ch_representative_transcript_gtf     // channel: [ val(meta), [ fai ] ]
    filtered_gtf                  = ch_filt_gtf                          // channel: [ val(meta), [ gtf ] ]
    seg_gtf                       = ch_seg_gtf                           // channel: [ val(meta), [ gtf ] ]
    regions_gtf                   = ch_regions_gtf                       // channel: [ val(meta), [ gtf ] ]
    regions_filt_gtf              = ch_regions_filt_gtf                  // channel: [ val(meta), [ gtf ] ]
    regions_resolved_gtf          = ch_regions_resolved_gtf              // channel: [ val(meta), [ gtf ] ]
    versions                      = ch_versions                          // channel: [ versions.yml ]
}
