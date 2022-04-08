package com.hartwig.hmftools.teal

import com.beust.jcommander.*
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.utils.FileWriterUtils

// only options that needed to be shown in validation are put here
private const val REF_SAMPLE = "-reference"
private const val TUMOR_SAMPLE = "-tumor"
private const val REF_WGS_METRICS = "-reference_wgs_metrics"
private const val TUMOR_WGS_METRICS = "-tumor_wgs_metrics"
private const val REF_BAM = "-reference_bam"
private const val TUMOR_BAM = "-tumor_bam"
private const val PURPLE = "-purple"
private const val COBALT = "-cobalt"
private const val REF_MEAN_READS_PER_KB = "-reference_mean_reads_per_kb"
private const val TUMOR_MEAN_READS_PER_KB = "-tumor_mean_reads_per_kb"

// NOTE: We use nullable string for required argument types, this is to avoid (default: <empty string>) getting
// printed in the usage by JCommander
data class TealCommonParams
    (
    @Parameter(names = [REF_SAMPLE], description = "ID of reference sample")
    var referenceSampleId: String? = null,

    @Parameter(names = [TUMOR_SAMPLE], description = "ID of tumor sample")
    var tumorSampleId: String? = null,

    @Parameter(names = [REF_BAM], description = "Path to reference bam/cram file")
    var referenceBamFile: String? = null,

    @Parameter(names = [TUMOR_BAM], description = "Path to tumor bam/cram file")
    var tumorBamFile: String? = null,

    @Parameter(names = ["-" + RefGenomeSource.REF_GENOME],
        description = "Path to reference genome fasta file if using CRAM files")
    var refGenomeFile: String? = null,

    @Parameter(names = ["-" + FileWriterUtils.OUTPUT_DIR],
        required = true,
        description = "Output directory")
    var outputDir: String? = null,

    @Parameter(names = ["-threads"], description = "Number of bam reader threads")
    var threadCount: Int = 1,

    @Parameter(names = ["-" + RefGenomeVersion.REF_GENOME_VERSION],
        description = RefGenomeVersion.REF_GENOME_VERSION_CFG_DESC,
        converter = RefGenomeVersionConverter::class)
    var refGenomeVersion: RefGenomeVersion = RefGenomeVersion.V37
    )
{
    fun getRefGenomeVersionStr() : String
    {
        return refGenomeVersion.identifier()
    }

    // we need to define a converter for ref genome version
    class RefGenomeVersionConverter : IStringConverter<RefGenomeVersion>
    {
        override fun convert(value: String): RefGenomeVersion
        {
            return RefGenomeVersion.from(value)
        }
    }

    fun referenceOnly() : Boolean
    {
        return referenceSampleId != null && tumorSampleId == null
    }

    fun tumorOnly() : Boolean
    {
        return referenceSampleId == null && tumorSampleId != null
    }

    // validate the params, throw ParameterException if fails
    @Throws(ParameterException::class)
    fun validate()
    {
        if (referenceSampleId == null && referenceBamFile != null)
        {
            throw ParameterException("$REF_SAMPLE is required when $REF_BAM(${referenceBamFile}) is specified")
        }

        if (referenceSampleId != null && referenceBamFile == null)
        {
            throw ParameterException("$REF_BAM is required when $REF_SAMPLE(${referenceSampleId}) is specified")
        }

        if (tumorSampleId == null && tumorBamFile != null)
        {
            throw ParameterException("$TUMOR_SAMPLE is required when ${TUMOR_BAM}(${tumorBamFile}) is specified")
        }

        if (tumorSampleId != null && tumorBamFile == null)
        {
            throw ParameterException("$TUMOR_BAM is required when $TUMOR_SAMPLE(${tumorSampleId}) is specified")
        }
    }
}

// NOTE: We use nullable string for required argument types, this is to avoid (default: <empty string>) getting
// printed in the usage by JCommander
data class TealPipelineParams
    (
    @Parameter(names = [PURPLE], description = "Path to PURPLE output directory")
    var purple: String? = null,

    @Parameter(names = [COBALT], required = true, description = "Path to COBALT output directory")
    var cobalt: String? = null,

    @Parameter(names = [REF_WGS_METRICS], description = "Path to reference WGS METRICS file")
    var referenceWgsMetrics: String? = null,

    @Parameter(names = [TUMOR_WGS_METRICS], description = "Path to tumor WGS METRICS file")
    var tumorWgsMetrics: String? = null,

    @ParametersDelegate
    val commonParams: TealCommonParams = TealCommonParams()
    )
{
    // validate the params, throw ParameterException if fails
    // NOTE we only validate optional arguments that could be missing
    @Throws(ParameterException::class)
    fun validate()
    {
        commonParams.validate()

        if (commonParams.referenceSampleId != null)
        {
            if (referenceWgsMetrics == null)
            {
                throw ParameterException("$REF_WGS_METRICS is required when $REF_SAMPLE(${commonParams.referenceSampleId}) is specified")
            }
        }

        if (commonParams.tumorSampleId != null)
        {
            if (purple == null)
            {
                throw ParameterException("$PURPLE is required when $TUMOR_SAMPLE(${commonParams.tumorSampleId}) is specified")
            }

            if (tumorWgsMetrics == null)
            {
                throw ParameterException("$TUMOR_WGS_METRICS is required when $TUMOR_SAMPLE(${commonParams.tumorSampleId}) is specified")
            }
        }
    }
}

// NOTE: We use nullable string for required argument types, this is to avoid (default: <empty string>) getting
// printed in the usage by JCommander
data class TealParams
    (
    @ParametersDelegate
    val commonParams: TealCommonParams = TealCommonParams(),

    @Parameter(names = ["-reference_duplicate_proportion"],
        description = "Proportion of reads that are marked duplicates in the reference sample BAM")
    var germlineDuplicateProportion: Double = 0.0,

    @Parameter(names = [REF_MEAN_READS_PER_KB], description = "Mean reads per KB of the reference sample")
    var germlineMeanReadsPerKb: Double? = null,

    @Parameter(names = ["-reference_gc50_reads_per_kb"],
        description = "GC 50 reads per KB of the reference sample. Defaults to mean reads per KB if not provided")
    var germlineGc50ReadsPerKb: Double? = null,

    @Parameter(names = ["-tumor_purity"],
        description = "Purity of the tumor sample")
    var tumorPurity: Double = 1.0,

    @Parameter(names = ["-tumor_ploidy"],
        description = "Ploidy of the tumor")
    var tumorPloidy: Double = 2.0,

    @Parameter(names = ["-tumor_duplicate_proportion"],
        description = "Proportion of reads that are marked duplicates in the tumor sample BAM")
    var tumorDuplicateProportion: Double = 0.0,

    @Parameter(names = [TUMOR_MEAN_READS_PER_KB], description = "Mean reads per KB of the tumor sample")
    var tumorMeanReadsPerKb: Double? = null,

    @Parameter(names = ["-tumor_gc50_reads_per_kb"], description = "GC 50 reads per KB. Defaults to mean reads per KB if not provided")
    var tumorGc50ReadsPerKb: Double? = null
)
{
    // validate the params, throw ParameterException if fails
    @Throws(ParameterException::class)
    fun validate()
    {
        commonParams.validate()

        if (commonParams.referenceSampleId != null)
        {
            if (germlineMeanReadsPerKb == null)
            {
                throw ParameterException("$REF_MEAN_READS_PER_KB is required when $REF_SAMPLE(${commonParams.referenceSampleId}) is specified")
            }
        }

        if (commonParams.tumorSampleId != null)
        {
            if (tumorMeanReadsPerKb == null)
            {
                throw ParameterException("$TUMOR_MEAN_READS_PER_KB is required when $TUMOR_SAMPLE(${commonParams.tumorSampleId}) is specified")
            }
        }
    }
}