package com.hartwig.hmftools.cider

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.perf.TaskExecutor
import com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE
import com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC
import com.hartwig.hmftools.common.utils.config.ConfigBuilder
import com.hartwig.hmftools.common.utils.file.FileWriterUtils
import org.apache.logging.log4j.LogManager


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
    val blast: String?,
    val blastDb: String?
)
{
    constructor(configBuilder: ConfigBuilder): this(
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
        blast = configBuilder.getValue(ARG_BLAST, null),
        blastDb = configBuilder.getValue(ARG_BLAST_DB, null)
    )
    {
        if ((blast == null) != (blastDb == null))
        {
            sLogger.error("invalid parameters: requires both -blast and -blast_db to be configBuilder.red together")
            throw IllegalArgumentException("Invalid blast configBuilder.ration")
        }
    }

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
        const val ARG_BLAST = "blast"
        const val ARG_BLAST_DB = "blast_db"

        private val sLogger = LogManager.getLogger(CiderParams::class.java)

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
            configBuilder.addPath(ARG_BLAST, false, "Location of blast installation")
            configBuilder.addPath(ARG_BLAST_DB, false, "Location of blast database")
        }
    }
}