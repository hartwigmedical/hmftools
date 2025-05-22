package com.hartwig.hmftools.teal

import com.beust.jcommander.*
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion
import com.hartwig.hmftools.common.perf.TaskExecutor
import com.hartwig.hmftools.common.utils.config.CommonConfig.*
import com.hartwig.hmftools.common.utils.config.ConfigBuilder
import com.hartwig.hmftools.common.utils.config.ConfigBuilder.getConfigDecimal
import com.hartwig.hmftools.common.utils.file.FileWriterUtils

private const val REF_WGS_METRICS = "reference_wgs_metrics"
private const val TUMOR_WGS_METRICS = "tumor_wgs_metrics"
private const val REF_MEAN_READ_DEPTH = "reference_mean_read_depth"
private const val TUMOR_MEAN_READ_DEPTH = "tumor_mean_read_depth"
private const val REFERENCE_DUPLICATE_PROPORTION = "reference_duplicate_proportion"
private const val REFERENCE_GC50_READ_DEPTH = "reference_gc50_read_depth"
private const val TUMOR_PURITY = "tumor_purity"
private const val TUMOR_PLOIDY = "tumor_ploidy"
private const val TUMOR_DUPLICATE_PROPORTION = "tumor_duplicate_proportion"
private const val TUMOR_GC50_READ_DEPTH = "tumor_gc50_read_depth"

// NOTE: We use nullable string for required argument types, this is to avoid (default: <empty string>) getting
// printed in the usage by JCommander
data class TealCommonParams(
    val referenceSampleId: String? = null,
    val tumorSampleId: String? = null,
    val referenceBamFile: String? = null,
    val tumorBamFile: String? = null,
    val refGenomeFile: String? = null,
    val outputDir: String? = null,
    val threadCount: Int = 1,
    val refGenomeVersion: RefGenomeVersion = RefGenomeVersion.V37
)
{
    constructor(configBuilder: ConfigBuilder): this(
        referenceSampleId = configBuilder.getValue(REFERENCE, null),
        tumorSampleId = configBuilder.getValue(TUMOR, null),
        referenceBamFile = configBuilder.getValue(REFERENCE_BAM, null),
        tumorBamFile = configBuilder.getValue(TUMOR_BAM, null),
        refGenomeFile = configBuilder.getValue(RefGenomeSource.REF_GENOME, null),
        outputDir = FileWriterUtils.parseOutputDir(configBuilder),
        threadCount = TaskExecutor.parseThreads(configBuilder),
        refGenomeVersion = RefGenomeVersion.from(configBuilder))

    fun getRefGenomeVersionStr() : String
    {
        return refGenomeVersion.identifier()
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
            throw ParameterException("$REFERENCE is required when $REFERENCE_BAM(${referenceBamFile}) is specified")
        }

        if (referenceSampleId != null && referenceBamFile == null)
        {
            throw ParameterException("$REFERENCE_BAM is required when $REFERENCE(${referenceSampleId}) is specified")
        }

        if (tumorSampleId == null && tumorBamFile != null)
        {
            throw ParameterException("$TUMOR is required when ${TUMOR_BAM}(${tumorBamFile}) is specified")
        }

        if (tumorSampleId != null && tumorBamFile == null)
        {
            throw ParameterException("$TUMOR_BAM is required when $TUMOR(${tumorSampleId}) is specified")
        }
    }
    
    // helper functions to standardise file names
    fun germlineTelbamPath() : String
    {
        return "${outputDir}/${referenceSampleId}.teal.telbam.bam"
    }

    fun tumorTelbamPath() : String
    {
        return "${outputDir}/${tumorSampleId}.teal.telbam.bam"
    }

    fun germlineTelLegnthTsvPath() : String
    {
        return "${outputDir}/${referenceSampleId}.teal.tellength.tsv"
    }

    fun tumorTelLengthTsvPath() : String
    {
        return "${outputDir}/${tumorSampleId}.teal.tellength.tsv"
    }

    fun breakEndTsvPath() : String
    {
        return "${outputDir}/${tumorSampleId}.teal.breakend.tsv.gz"
    }

    companion object
    {
        fun registerConfig(configBuilder: ConfigBuilder)
        {
            configBuilder.addConfigItem(REFERENCE, REFERENCE_DESC)
            configBuilder.addPath(REFERENCE_BAM, false, REFERENCE_BAM_DESC)
            configBuilder.addConfigItem(TUMOR, TUMOR_DESC)
            configBuilder.addPath(TUMOR_BAM, false, TUMOR_BAM_DESC)
            addRefGenomeFile(configBuilder, false)
            FileWriterUtils.addOutputDir(configBuilder)
            TaskExecutor.addThreadOptions(configBuilder)
            RefGenomeSource.addRefGenomeVersion(configBuilder)
        }
    }
}

// NOTE: We use nullable string for required argument types, this is to avoid (default: <empty string>) getting
// printed in the usage by JCommander
data class TealPipelineParams(
    val purple: String? = null,
    val cobalt: String? = null,
    val referenceWgsMetrics: String? = null,
    val tumorWgsMetrics: String? = null,
    val commonParams: TealCommonParams = TealCommonParams()
    )
{
    constructor(configBuilder: ConfigBuilder): this(
        purple = configBuilder.getValue("purple"),
        cobalt = configBuilder.getValue("cobalt"),
        referenceWgsMetrics = configBuilder.getValue(REF_WGS_METRICS, null),
        tumorWgsMetrics = configBuilder.getValue(TUMOR_WGS_METRICS, null),
        commonParams = TealCommonParams(configBuilder))

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
                throw ParameterException("$REF_WGS_METRICS is required when $REFERENCE(${commonParams.referenceSampleId}) is specified")
            }
        }

        if (commonParams.tumorSampleId != null)
        {
            if (purple == null)
            {
                throw ParameterException("$PURPLE_DIR_CFG is required when $TUMOR(${commonParams.tumorSampleId}) is specified")
            }

            if (tumorWgsMetrics == null)
            {
                throw ParameterException("$TUMOR_WGS_METRICS is required when $TUMOR(${commonParams.tumorSampleId}) is specified")
            }
        }
    }

    companion object
    {
        fun registerConfig(configBuilder: ConfigBuilder)
        {
            configBuilder.addPath("purple", false, PURPLE_DIR_DESC)
            configBuilder.addPath("cobalt", true, COBALT_DIR_CFG)
            configBuilder.addPath(REF_WGS_METRICS, false, "Path to reference BAM METRICS file")
            configBuilder.addPath(TUMOR_WGS_METRICS, false, "Path to tumor BAM METRICS file")
            TealCommonParams.registerConfig(configBuilder)
        }
    }
}

// NOTE: We use nullable string for required argument types, this is to avoid (default: <empty string>) getting
// printed in the usage by JCommander
data class TealParams(
    val commonParams: TealCommonParams = TealCommonParams(),
    var germlineDuplicateProportion: Double = 0.0,
    var germlineMeanReadDepth: Double? = null,
    var germlineGc50ReadDepth: Double? = null,
    var tumorPurity: Double = 1.0,
    var tumorPloidy: Double = 2.0,
    var tumorDuplicateProportion: Double = 0.0,
    var tumorMeanReadDepth: Double? = null,
    var tumorGc50ReadDepth: Double? = null
)
{
    constructor(configBuilder: ConfigBuilder): this(
        commonParams = TealCommonParams(configBuilder),
        germlineDuplicateProportion = getConfigDecimal(configBuilder, REFERENCE_DUPLICATE_PROPORTION, 0.0),
        germlineMeanReadDepth = if (configBuilder.hasValue(REF_MEAN_READ_DEPTH))
            configBuilder.getDecimal(REF_MEAN_READ_DEPTH) else null,
        germlineGc50ReadDepth = if (configBuilder.hasValue(REFERENCE_GC50_READ_DEPTH))
            configBuilder.getDecimal(REFERENCE_GC50_READ_DEPTH) else null,
        tumorPurity = configBuilder.getDecimal(TUMOR_PURITY),
        tumorPloidy = configBuilder.getDecimal(TUMOR_PLOIDY),
        tumorDuplicateProportion = getConfigDecimal(configBuilder, TUMOR_DUPLICATE_PROPORTION, 0.0),
        tumorMeanReadDepth = if (configBuilder.hasValue(TUMOR_MEAN_READ_DEPTH))
            configBuilder.getDecimal(TUMOR_MEAN_READ_DEPTH) else null,
        tumorGc50ReadDepth = if (configBuilder.hasValue(TUMOR_GC50_READ_DEPTH))
            configBuilder.getDecimal(TUMOR_GC50_READ_DEPTH) else null,
    )

    // validate the params, throw ParameterException if fails
    @Throws(ParameterException::class)
    fun validate()
    {
        commonParams.validate()

        if (commonParams.referenceSampleId != null)
        {
            if (germlineMeanReadDepth == null)
            {
                throw ParameterException("$REF_MEAN_READ_DEPTH is required when $REFERENCE(${commonParams.referenceSampleId}) is specified")
            }
        }

        if (commonParams.tumorSampleId != null)
        {
            if (tumorMeanReadDepth == null)
            {
                throw ParameterException("$TUMOR_MEAN_READ_DEPTH is required when $TUMOR(${commonParams.tumorSampleId}) is specified")
            }
        }
    }

    companion object
    {
        fun registerConfig(configBuilder: ConfigBuilder)
        {
            TealCommonParams.registerConfig(configBuilder)
            configBuilder.addConfigItem(
                REFERENCE_DUPLICATE_PROPORTION,
                "Proportion of reads that are marked duplicates in the reference sample BAM")
            configBuilder.addConfigItem(REF_MEAN_READ_DEPTH, "Mean read depth of the reference sample")
            configBuilder.addConfigItem(REFERENCE_GC50_READ_DEPTH, "GC 50 read depth of the reference sample")
            configBuilder.addDecimal(TUMOR_PURITY, "Purity of the tumor sample", 1.0)
            configBuilder.addDecimal(TUMOR_PLOIDY, "Ploidy of the tumor", 2.0)
            configBuilder.addDecimal(TUMOR_DUPLICATE_PROPORTION, "Proportion of reads that are marked duplicates in the tumor sample BAM", 0.0)
            configBuilder.addConfigItem(TUMOR_MEAN_READ_DEPTH, "Mean read depth of the tumor sample")
            configBuilder.addConfigItem(TUMOR_GC50_READ_DEPTH, "GC 50 read depth of the tumor sample")
        }
    }
}