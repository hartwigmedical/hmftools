package com.hartwig.hmftools.cider

import com.hartwig.hmftools.common.bwa.BwaUtils
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.perf.TaskExecutor
import com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE
import com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC
import com.hartwig.hmftools.common.utils.config.ConfigBuilder
import com.hartwig.hmftools.common.utils.file.FileWriterUtils


data class CiderParams(
    val sampleId: String,
    val bamPath: String,
    val refGenomePath: String?,
    val outputDir: String,
    val threadCount: Int,
    val refGenomeVersion: RefGenomeVersion,
    val approxMaxFragmentLength: Int,
    val minBaseQuality: Int,
    val writeFilteredBam: Boolean,
    val reportMatchRefSeq: Boolean,
    val numBasesToTrim: Int,
    val maxLowQualBaseFraction: Double,
    val maxReadCountPerGene: Int,
    val primerCsv: String?,
    val primerMismatchMax: Int,
    val bwaLibPath: String?,
    val alignmentRefGenomePath: String?,
    val alignmentBwaIndexImagePath: String?
)
{
    companion object
    {
        const val DEFAULT_MAX_FRAGMENT_LENGTH = 1000
        const val DEFAULT_MAX_READ_COUNT_PER_GENE = 600_000

        // maximum proportion of read bases that are low quality
        const val DEFAULT_MAX_LOW_QUAL_BASES_FRACTION: Double = 0.1
        const val DEFAULT_MIN_BASE_QUALITY = 25

        const val ARG_BAM = "bam"
        const val ARG_MAX_FRAGMENT_LENGTH = "max_fragment_length"
        const val ARG_MIN_BASE_QUALITY = "min_base_quality"
        const val ARG_WRITE_CIDER_BAM = "write_cider_bam"
        const val ARG_REPORT_MATCH_REF_SEQ = "report_match_ref_seq"
        const val ARG_NUM_TRIM_BASES = "num_trim_bases"
        const val ARG_MAX_LOW_QUAL_BASE_FRACTION = "max_low_qual_base_fraction"
        const val ARG_MAX_READS_PER_GENE = "max_reads_per_gene"
        const val ARG_PRIMER_CSV = "primer_csv"
        const val ARG_PRIMER_MISMATCH_MAX = "primer_mismatch_max"
        const val ARG_ALIGNMENT_REF_GENOME = "alignment_ref_genome"
        const val ARG_BWA_INDEX_IMAGE_FILE = "bwa_index_image"

        fun fromConfigBuilder(configBuilder: ConfigBuilder): CiderParams {
            val alignmentRefGenome = configBuilder.getValue(ARG_ALIGNMENT_REF_GENOME)
            val defaultAlignmentBwaIndex = if (alignmentRefGenome == null) { null } else { "$alignmentRefGenome.img" }
            return CiderParams(
                sampleId = configBuilder.getValue(SAMPLE),
                bamPath = configBuilder.getValue(ARG_BAM),
                refGenomePath = configBuilder.getValue(RefGenomeSource.REF_GENOME, null),
                outputDir = FileWriterUtils.parseOutputDir(configBuilder),
                threadCount = TaskExecutor.parseThreads(configBuilder),
                approxMaxFragmentLength = configBuilder.getInteger(ARG_MAX_FRAGMENT_LENGTH),
                refGenomeVersion = RefGenomeVersion.from(configBuilder),
                minBaseQuality = configBuilder.getInteger(ARG_MIN_BASE_QUALITY),
                writeFilteredBam = configBuilder.hasFlag(ARG_WRITE_CIDER_BAM),
                reportMatchRefSeq = configBuilder.hasFlag(ARG_REPORT_MATCH_REF_SEQ),
                numBasesToTrim = configBuilder.getInteger(ARG_NUM_TRIM_BASES),
                maxLowQualBaseFraction = configBuilder.getDecimal(ARG_MAX_LOW_QUAL_BASE_FRACTION),
                maxReadCountPerGene = configBuilder.getInteger(ARG_MAX_READS_PER_GENE),
                primerCsv = configBuilder.getValue(ARG_PRIMER_CSV),
                primerMismatchMax = configBuilder.getInteger(ARG_PRIMER_MISMATCH_MAX),
                bwaLibPath = configBuilder.getValue(BwaUtils.BWA_LIB_PATH),
                alignmentRefGenomePath = alignmentRefGenome,
                alignmentBwaIndexImagePath = configBuilder.getValue(ARG_BWA_INDEX_IMAGE_FILE, defaultAlignmentBwaIndex)
            )
        }

        fun registerConfig(configBuilder: ConfigBuilder)
        {
            configBuilder.addConfigItem(SAMPLE, SAMPLE_DESC)
            configBuilder.addPath(ARG_BAM, true, "input BAM / CRAM file")
            addRefGenomeFile(configBuilder, false)
            FileWriterUtils.addOutputDir(configBuilder)
            TaskExecutor.addThreadOptions(configBuilder)
            RefGenomeSource.addRefGenomeVersion(configBuilder)
            configBuilder.addInteger(ARG_MAX_FRAGMENT_LENGTH, "Approximate maximum fragment length", DEFAULT_MAX_FRAGMENT_LENGTH)
            configBuilder.addInteger(ARG_MIN_BASE_QUALITY, "Minimum quality for a base to be considered", DEFAULT_MIN_BASE_QUALITY)
            configBuilder.addFlag(ARG_WRITE_CIDER_BAM, "Write a output BAM file containing all CDR3 reads")
            configBuilder.addFlag(ARG_REPORT_MATCH_REF_SEQ, "Report VDJ sequences that match reference genome")
            configBuilder.addInteger(ARG_NUM_TRIM_BASES, "Number of bases to trim on each side of reads", 0)
            configBuilder.addDecimal(
                ARG_MAX_LOW_QUAL_BASE_FRACTION,
                "Maximum fraction of bases in a read that can be low quality. Reads that exceed this limit are discarded",
                DEFAULT_MAX_LOW_QUAL_BASES_FRACTION
            )
            configBuilder.addInteger(
                ARG_MAX_READS_PER_GENE,
                "Maximum number of reads per gene. If number of reads exceed this limit, they are downsampled.",
                DEFAULT_MAX_READ_COUNT_PER_GENE
            )
            configBuilder.addPath(ARG_PRIMER_CSV, false, "Path to csv file containing primers")
            configBuilder.addInteger(
                ARG_PRIMER_MISMATCH_MAX,
                "Maximum number of mismatch bases for matching primer sequence",
                0
            )
            configBuilder.addPath(BwaUtils.BWA_LIB_PATH, false, BwaUtils.BWA_LIB_PATH_DESC)
            configBuilder.addPath(ARG_ALIGNMENT_REF_GENOME, false, "Reference genome FASTA for alignment")
            configBuilder.addPath(ARG_BWA_INDEX_IMAGE_FILE, false, "Reference genome BWA-MEM index GATK image file")
        }
    }
}