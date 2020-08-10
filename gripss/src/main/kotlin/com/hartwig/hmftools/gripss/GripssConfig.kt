package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.common.cli.Configs.defaultDoubleValue
import com.hartwig.hmftools.common.cli.Configs.defaultIntValue
import htsjdk.samtools.reference.IndexedFastaSequenceFile
import htsjdk.samtools.util.Interval
import htsjdk.samtools.util.Locatable
import org.apache.commons.cli.CommandLine
import org.apache.commons.cli.Option
import org.apache.commons.cli.Options
import org.apache.commons.cli.ParseException
import org.apache.logging.log4j.util.Strings
import java.io.File
import java.io.IOException

const val INPUT_VCF_OPTION = "input_vcf"
const val OUTPUT_VCF_OPTION = "output_vcf"
const val REF_GENOME_OPTION = "ref_genome"
const val SINGLE_PON_OPTION = "breakend_pon"
const val PAIRED_PON_OPTION = "breakpoint_pon"
const val PAIRED_HOTSPOT_OPTION = "breakpoint_hotspot"
const val REFERENCE = "reference"
const val TUMOR = "tumor"

data class GripssConfig(
        val inputVcf: String,
        val outputVcf: String,
        val singlePonFile: String,
        val pairedPonFile: String,
        val pairedHotspotFile: String,
        val refGenome: String,
        val reference: String,
        val tumor: String,
        val filterConfig: GripssFilterConfig) {

    companion object {
        fun createOptions(): Options {
            val options = Options()
            options.addOption(requiredOption(INPUT_VCF_OPTION, "Path to GRIDSS VCF input"))
            options.addOption(requiredOption(OUTPUT_VCF_OPTION, "Path to output VCF"))
            options.addOption(requiredOption(SINGLE_PON_OPTION, "Single breakend pon bed file"))
            options.addOption(requiredOption(PAIRED_PON_OPTION, "Paired breakpoint pon bedpe file"))
            options.addOption(requiredOption(PAIRED_HOTSPOT_OPTION, "Paired breakpoint hotspot bedpe file"))
            options.addOption(requiredOption(REF_GENOME_OPTION, "Ref genome"))
            options.addOption(Option(REFERENCE, true, "Optional name of reference sample [first sample in VCF]"))
            options.addOption(Option(TUMOR, true, "Optional name of tumor sample [second sample in VCF]"))
            GripssFilterConfig.createOptions().options.forEach { options.addOption(it) }
            return options
        }

        @Throws(ParseException::class, IOException::class)
        fun createConfig(cmd: CommandLine): GripssConfig {
            val inputVcf = requiredFile(cmd, INPUT_VCF_OPTION)
            val refGenome = requiredFile(cmd, REF_GENOME_OPTION)
            val singlePon = requiredFile(cmd, SINGLE_PON_OPTION)
            val pairedPon = requiredFile(cmd, PAIRED_PON_OPTION)
            val pairedHotspot = requiredFile(cmd, PAIRED_HOTSPOT_OPTION)
            val outputVcf = cmd.getOptionValue(OUTPUT_VCF_OPTION)
            val reference = cmd.getOptionValue(REFERENCE, Strings.EMPTY)
            val tumor = cmd.getOptionValue(TUMOR, Strings.EMPTY)

            val outputDir = File(outputVcf).absoluteFile.parentFile
            if (!outputDir.exists() && !outputDir.mkdirs()) {
                throw IOException("Unable to write to directory ${outputDir.absolutePath}")
            }

            val isGRCh38 = isGRCh38(refGenome)
            val filterConfig = GripssFilterConfig.createConfig(cmd, isGRCh38)

            return GripssConfig(inputVcf, outputVcf, singlePon, pairedPon, pairedHotspot, refGenome, reference, tumor, filterConfig)
        }

        private fun isGRCh38(refGenome: String): Boolean {
            IndexedFastaSequenceFile(File(refGenome)).use {
                val dict = it.sequenceDictionary
                return dict.getSequence("chr2")?.sequenceLength == 242_193_529
            }
        }

        @Throws(IOException::class)
        private fun requiredFile(cmd: CommandLine, argument: String): String {
            val result = cmd.getOptionValue(argument)
            if (!File(result).exists()) {
                throw IOException("Unable to read file $result")
            }

            return result
        }

        private fun requiredOption(argument: String, description: String): Option {
            val result = Option(argument, true, description)
            result.isRequired = true
            return result
        }
    }
}

private const val HARD_MIN_TUMOR_QUAL_OPTION = "hard_min_tumor_qual"
private const val HARD_MAX_NORMAL_ABSOLUTE_SUPPORT_OPTION = "hard_max_normal_absolute_support"
private const val HARD_MAX_NORMAL_RELATIVE_SUPPORT_OPTION = "hard_max_normal_relative_support"
private const val SOFT_MAX_NORMAL_RELATIVE_SUPPORT_OPTION = "soft_max_normal_relative_support"
private const val MIN_NORMAL_COVERAGE_OPTION = "min_normal_coverage"
private const val MIN_TUMOR_AF_OPTION = "min_tumor_af"
private const val MAX_SHORT_STRAND_BIAS_OPTION = "max_short_strand_bias"
private const val MIN_QUAL_BREAK_END_OPTION = "min_qual_break_end"
private const val MIN_QUAL_BREAK_POINT_OPTION = "min_qual_break_point"
private const val MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION = "min_qual_rescue_mobile_element_insertion"
private const val MAX_HOM_LENGTH_SHORT_INV_OPTION = "max_hom_length_short_inv"
private const val MAX_INEXACT_HOM_LENGTH_SHORT_DEL_OPTION = "max_inexact_hom_length_short_del"
private const val MIN_LENGTH_OPTION = "min_length"

data class GripssFilterConfig(
        val hardMinTumorQual: Int,
        val hardMaxNormalAbsoluteSupport: Int,
        val hardMaxNormalRelativeSupport: Double,
        val softMaxNormalRelativeSupport: Double,
        val minNormalCoverage: Int,
        val minTumorAF: Double,
        val maxShortStrandBias: Double,
        val minQualBreakEnd: Int,
        val minQualBreakPoint: Int,
        val minQualRescueMobileElementInsertion: Int,
        val maxHomLengthShortInversion: Int,
        val maxInexactHomLengthShortDel: Int,
        val minLength: Int,
        val polyGCRegion: Locatable) {

    companion object {

        private val linc00486Definition37 = Interval("2", 33_141_260, 33_141_700)
        private val linc00486Definition38 = Interval("2", 32_916_190, 32_916_630)

        fun createOptions(): Options {
            val defaultConfig = default()
            val options = Options()
            options.addOption(HARD_MIN_TUMOR_QUAL_OPTION, "Hard min tumor qual [${defaultConfig.hardMinTumorQual}]")
            options.addOption(HARD_MAX_NORMAL_ABSOLUTE_SUPPORT_OPTION, "Hard max normal absolute support [${defaultConfig.hardMaxNormalAbsoluteSupport}]")
            options.addOption(HARD_MAX_NORMAL_RELATIVE_SUPPORT_OPTION, "Hard max normal relative support [${defaultConfig.hardMaxNormalRelativeSupport}]")
            options.addOption(SOFT_MAX_NORMAL_RELATIVE_SUPPORT_OPTION, "Soft max normal relative support [${defaultConfig.softMaxNormalRelativeSupport}]")
            options.addOption(MIN_NORMAL_COVERAGE_OPTION, "Min normal coverage [${defaultConfig.minNormalCoverage}]")
            options.addOption(MIN_TUMOR_AF_OPTION, "Min tumor allelic frequency [${defaultConfig.minTumorAF}]")
            options.addOption(MAX_SHORT_STRAND_BIAS_OPTION, "Max short strand bias [${defaultConfig.maxShortStrandBias}]")
            options.addOption(MIN_QUAL_BREAK_END_OPTION, "Min qual break end [${defaultConfig.minQualBreakEnd}]")
            options.addOption(MIN_QUAL_BREAK_POINT_OPTION, "Min qual break point [${defaultConfig.minQualBreakPoint}]")
            options.addOption(MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION, "Min qual rescue mobile element insertions [${defaultConfig.minQualRescueMobileElementInsertion}]")
            options.addOption(MAX_HOM_LENGTH_SHORT_INV_OPTION, "Max homology length short inversion [${defaultConfig.maxHomLengthShortInversion}]")
            options.addOption(MAX_INEXACT_HOM_LENGTH_SHORT_DEL_OPTION, "Max inexact homology length short del [${defaultConfig.maxInexactHomLengthShortDel}]")
            options.addOption(MIN_LENGTH_OPTION, "Min length [${defaultConfig.minLength}]")

            return options
        }

        fun createConfig(cmd: CommandLine, isGRCh38: Boolean): GripssFilterConfig {
            val defaultConfig = default()
            val minTumorQual = defaultIntValue(cmd, HARD_MIN_TUMOR_QUAL_OPTION, defaultConfig.hardMinTumorQual)
            val maxNormalAbsoluteSupport = defaultIntValue(cmd, HARD_MAX_NORMAL_ABSOLUTE_SUPPORT_OPTION, defaultConfig.hardMaxNormalAbsoluteSupport)
            val hardMaxNormalRelativeSupport = defaultDoubleValue(cmd, HARD_MAX_NORMAL_RELATIVE_SUPPORT_OPTION, defaultConfig.hardMaxNormalRelativeSupport)
            val softMaxNormalRelativeSupport = defaultDoubleValue(cmd, SOFT_MAX_NORMAL_RELATIVE_SUPPORT_OPTION, defaultConfig.softMaxNormalRelativeSupport)
            val minNormalCoverage = defaultIntValue(cmd, MIN_NORMAL_COVERAGE_OPTION, defaultConfig.minNormalCoverage)
            val minTumorAF = defaultDoubleValue(cmd, MIN_TUMOR_AF_OPTION, defaultConfig.minTumorAF)
            val maxShortStrandBias = defaultDoubleValue(cmd, MAX_SHORT_STRAND_BIAS_OPTION, defaultConfig.maxShortStrandBias)
            val minQualBreakEnd = defaultIntValue(cmd, MIN_QUAL_BREAK_END_OPTION, defaultConfig.minQualBreakEnd)
            val minQualBreakPoint = defaultIntValue(cmd, MIN_QUAL_BREAK_POINT_OPTION, defaultConfig.minQualBreakPoint)
            val minQualRescueMobileElementInserions = defaultIntValue(cmd, MIN_QUAL_RESCUE_MOBILE_ELEMENT_INSERTION, defaultConfig.minQualRescueMobileElementInsertion)
            val maxHomLengthShortInversion = defaultIntValue(cmd, MAX_HOM_LENGTH_SHORT_INV_OPTION, defaultConfig.maxHomLengthShortInversion)
            val maxInexactHomLengthShortDel = defaultIntValue(cmd, MAX_INEXACT_HOM_LENGTH_SHORT_DEL_OPTION, defaultConfig.maxInexactHomLengthShortDel)
            val minLength = defaultIntValue(cmd, MIN_LENGTH_OPTION, defaultConfig.minLength)
            val linc00486Definition = if (isGRCh38) linc00486Definition38 else linc00486Definition37

            return GripssFilterConfig(
                    minTumorQual,
                    maxNormalAbsoluteSupport,
                    hardMaxNormalRelativeSupport,
                    softMaxNormalRelativeSupport,
                    minNormalCoverage,
                    minTumorAF,
                    maxShortStrandBias,
                    minQualBreakEnd,
                    minQualBreakPoint,
                    minQualRescueMobileElementInserions,
                    maxHomLengthShortInversion,
                    maxInexactHomLengthShortDel,
                    minLength,
                    linc00486Definition)
        }

        fun default(): GripssFilterConfig {
            return GripssFilterConfig(
                    100,
                    3,
                    0.06,
                    0.03,
                    8,
                    0.005,
                    0.95,
                    1000,
                    400,
                    1000,
                    6,
                    5,
                    32,
                    linc00486Definition37)
        }
    }
}