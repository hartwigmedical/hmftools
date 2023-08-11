package com.hartwig.hmftools.cider

import com.beust.jcommander.Parameter
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.utils.config.RefGenomeVersionConverter
import org.apache.logging.log4j.LogManager

class CiderParams
{
    @Parameter(names = ["-sample"], description = "Name of sample")
    lateinit var sampleId: String

    @Parameter(names = ["-bam"], description = "Path to indexed bam/cram file")
    lateinit var bamPath: String

    @Parameter(
        names = ["-" + RefGenomeSource.REF_GENOME],
        description = "Path to the reference genome fasta file. Required only when using CRAM files."
    )
    var refGenomePath: String? = null

    @Parameter(
        names = ["-output_dir"],
        required = true,
        description = "Path to the output directory. "
                + "This directory will be created if it does not already exist."
    )
    lateinit var outputDir: String

    @Parameter(names = ["-threads"], description = "Number of threads")
    var threadCount = DEFAULT_THREADS

    @Parameter(names = ["-max_fragment_length"], description = "Approximate maximum fragment length")
    var approxMaxFragmentLength = DEFAULT_MAX_FRAGMENT_LENGTH

    @Parameter(
        names = ["-" + RefGenomeVersion.REF_GENOME_VERSION],
        required = true,
        description = RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC,
        converter = RefGenomeVersionConverter::class
    )
    lateinit var refGenomeVersion: RefGenomeVersion

    @Parameter(names = ["-min_base_quality"], description = "Minimum quality for a base to be considered")
    var minBaseQuality = 25

    @Parameter(names = ["-write_cider_bam"], description = "Write a output BAM file containing all CDR3 reads")
    var writeFilteredBam = false

    @Parameter(names = ["-report_match_ref_seq"], description = "Report VDJ sequences that match reference genome")
    var reportMatchRefSeq = false

    @Parameter(names = ["-num_trim_bases"], description = "Number of bases to trim on each side of reads")
    var numBasesToTrim = 0

    @Parameter(names = ["-max_low_qual_base_fraction"],
        description = "Maximum fraction of bases in a read that can be low quality." +
                "Reads that exceed this limit are discarded")
    var maxLowQualBaseFraction = MAX_LOW_QUAL_BASES_FRACTION

    @Parameter(names = ["-max_reads_per_gene"],
        description = "Maximum number of reads per gene. If number of reads exceed this limit, they are downsampled.")
    var maxReadCountPerGene = DEFAULT_MAX_READ_COUNT_PER_GENE

    @Parameter(names = ["-primer_csv"], description = "Path to csv file containing primers")
    var primerCsv: String? = null

    @Parameter(names = ["-primer_mismatch_max"], description = "Maximum number of mismatch bases for matching primer sequence")
    var primerMismatchMax: Int = 0

    @Parameter(names = ["-blast"], description = "Location of blast installation")
    var blast: String? = null

    @Parameter(names = ["-blast_db"], description = "Location of blast database")
    var blastDb: String? = null

    val isValid: Boolean get()
    {
        if (blast != null && blastDb == null)
        {
            sLogger.error("invalid parameters: requires -blast_db if -blast is configured")
            return false
        }
        if (blast == null && blastDb != null)
        {
            sLogger.error("invalid parameters: requires -blast if -blast_db is configured")
            return false
        }
        return true
    }

    companion object
    {
        const val DEFAULT_THREADS = 1
        const val DEFAULT_MAX_FRAGMENT_LENGTH = 1000
        const val DEFAULT_MAX_READ_COUNT_PER_GENE = 600_000
        // maximum proportion of read bases that are low quality
        const val MAX_LOW_QUAL_BASES_FRACTION: Double = 0.1

        val sLogger = LogManager.getLogger(CiderParams::class.java)
    }
}