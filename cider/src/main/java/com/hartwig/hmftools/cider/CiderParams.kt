package com.hartwig.hmftools.cider

import com.beust.jcommander.Parameter
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.utils.config.RefGenomeVersionConverter
import htsjdk.samtools.ValidationStringency

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

    @Parameter(names = ["-validation_stringency"], description = "SAM validation strategy")
    var stringency = ValidationStringency.DEFAULT_STRINGENCY

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

    @Parameter(names = ["-" + EnsemblDataCache.ENSEMBL_DATA_DIR],
               required = true,
               description = EnsemblDataCache.ENSEMBL_DATA_DIR_CFG)
    lateinit var ensemblDataDir: String

    @Parameter(names = ["-min_base_quality"], description = "Minimum quality for a base to be considered")
    var minBaseQuality = 25

    @Parameter(names = ["-write_cider_bam"], description = "Write a output BAM file containing all CDR3 reads")
    var writeFilteredBam = false

    @Parameter(names = ["-num_trim_bases"], description = "Number of bases to trim on each side of reads")
    var numBasesToTrim = 0

    @Parameter(names = ["-primer_csv"], description = "Path to csv file containing primers")
    var primerCsv: String? = null

    @Parameter(names = ["-primer_mismatch_max"], description = "Maximum number of mismatch bases for matching primer sequence")
    var primerMismatchMax: Int = 0

    val isValid: Boolean get() = true

    companion object
    {
        const val DEFAULT_THREADS = 1
        const val DEFAULT_MAX_FRAGMENT_LENGTH = 1000
    }
}