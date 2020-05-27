package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.common.cli.Configs.defaultDoubleValue
import com.hartwig.hmftools.common.cli.Configs.defaultIntValue
import org.apache.commons.cli.CommandLine
import org.apache.commons.cli.Option
import org.apache.commons.cli.Options
import org.apache.commons.cli.ParseException
import java.io.File
import java.io.IOException

const val INPUT_VCF_OPTION = "input_vcf"
const val OUTPUT_VCF_OPTION = "output_vcf"
const val REF_GENOME_OPTION = "ref_genome"
const val SINGLE_PON_OPTION = "breakend_pon"
const val PAIRED_PON_OPTION = "breakpoint_pon"
const val PAIRED_HOTSPOT_OPTION = "breakpoint_hotspot"

data class GripssConfig(
        val inputVcf: String,
        val outputVcf: String,
        val singlePonFile: String,
        val pairedPonFile: String,
        val pairedHotspotFile: String,
        val refGenome: String,
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
            GripssFilterConfig.createOptions().options.forEach { options.addOption(it) }
            return options
        }

        @Throws(ParseException::class, IOException::class)
        fun createConfig(cmd: CommandLine): GripssConfig {
            val filterConfig = GripssFilterConfig.createConfig(cmd)
            val inputVcf = requiredFile(cmd, INPUT_VCF_OPTION)
            val refGenome = requiredFile(cmd, REF_GENOME_OPTION)
            val singlePon = requiredFile(cmd, SINGLE_PON_OPTION)
            val pairedPon = requiredFile(cmd, PAIRED_PON_OPTION)
            val pairedHotspot = requiredFile(cmd, PAIRED_HOTSPOT_OPTION)
            val outputVcf = cmd.getOptionValue(OUTPUT_VCF_OPTION)

            val outputDir = File(outputVcf).parentFile
            if (!outputDir.exists() && !outputDir.mkdirs()) {
                throw IOException("Unable to write to directory ${outputDir.absolutePath}")
            }

            return GripssConfig(inputVcf, outputVcf, singlePon, pairedPon, pairedHotspot, refGenome, filterConfig)
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

private const val MAX_NORMAL_SUPPORT_PROPORTION_OPTION = "max_normal_support_proportion"
private const val MIN_NORMAL_COVERAGE_OPTION = "min_normal_coverage"
private const val MIN_TUMOR_AF_OPTION = "min_tumor_af"
private const val MAX_SHORT_STRAND_BIAS_OPTION = "max_short_strand_bias"
private const val MIN_QUAL_BREAK_END_OPTION = "min_qual_break_end"
private const val MIN_QUAL_BREAK_POINT_OPTION = "min_qual_break_point"
private const val MAX_HOM_LENGTH_SHORT_INV_OPTION = "max_hom_length_short_inv"
private const val MAX_INEXACT_HOM_LENGTH_OPTION = "max_inexact_hom_length"
private const val MAX_INEXACT_HOM_LENGTH_SHORT_DEL_OPTION = "max_inexact_hom_length_short_del"
private const val MIN_LENGTH_OPTION = "min_length"

data class GripssFilterConfig(
        val maxNormalSupportProportion: Double,
        val minNormalCoverage: Int,
        val minTumorAF: Double,
        val maxShortStrandBias: Double,
        val minQualBreakEnd: Int,
        val minQualBreakPoint: Int,
        val maxHomLengthShortInversion: Int,
        val maxInexactHomLength: Int,
        val maxInexactHomLengthShortDel: Int,
        val minLength: Int) {

    companion object {
        fun createOptions(): Options {
            val defaultConfig = default()
            val options = Options()
            options.addOption(MAX_NORMAL_SUPPORT_PROPORTION_OPTION, "Max normal support [${defaultConfig.maxNormalSupportProportion}]")
            options.addOption(MIN_NORMAL_COVERAGE_OPTION, "Min normal coverage [${defaultConfig.minNormalCoverage}]")
            options.addOption(MIN_TUMOR_AF_OPTION, "Min tumor allelic frequency [${defaultConfig.minTumorAF}]")
            options.addOption(MAX_SHORT_STRAND_BIAS_OPTION, "Max short strand bias [${defaultConfig.maxShortStrandBias}]")
            options.addOption(MIN_QUAL_BREAK_END_OPTION, "Min qual break end [${defaultConfig.minQualBreakEnd}]")
            options.addOption(MIN_QUAL_BREAK_POINT_OPTION, "Min qual break point [${defaultConfig.minQualBreakPoint}]")
            options.addOption(MAX_HOM_LENGTH_SHORT_INV_OPTION, "Max homology length short inversion [${defaultConfig.maxHomLengthShortInversion}]")
            options.addOption(MAX_INEXACT_HOM_LENGTH_OPTION, "Max inexact homology length [${defaultConfig.maxInexactHomLength}]")
            options.addOption(MAX_INEXACT_HOM_LENGTH_SHORT_DEL_OPTION, "Max inexact homology length short del [${defaultConfig.maxInexactHomLengthShortDel}]")
            options.addOption(MIN_LENGTH_OPTION, "Min length [${defaultConfig.minLength}]")

            return options
        }

        fun createConfig(cmd: CommandLine): GripssFilterConfig {
            val defaultConfig = default()
            val maxNormalSupportProportion = defaultDoubleValue(cmd, MAX_NORMAL_SUPPORT_PROPORTION_OPTION, defaultConfig.maxNormalSupportProportion)
            val minNormalCoverage = defaultIntValue(cmd, MIN_NORMAL_COVERAGE_OPTION, defaultConfig.minNormalCoverage)
            val minTumorAF = defaultDoubleValue(cmd, MIN_TUMOR_AF_OPTION, defaultConfig.minTumorAF)
            val maxShortStrandBias = defaultDoubleValue(cmd, MAX_SHORT_STRAND_BIAS_OPTION, defaultConfig.maxShortStrandBias)
            val minQualBreakEnd = defaultIntValue(cmd, MIN_QUAL_BREAK_END_OPTION, defaultConfig.minQualBreakEnd)
            val minQualBreakPoint = defaultIntValue(cmd, MIN_QUAL_BREAK_POINT_OPTION, defaultConfig.minQualBreakPoint)
            val maxHomLengthShortInversion = defaultIntValue(cmd, MAX_HOM_LENGTH_SHORT_INV_OPTION, defaultConfig.maxHomLengthShortInversion)
            val maxInexactHomLength = defaultIntValue(cmd, MAX_INEXACT_HOM_LENGTH_OPTION, defaultConfig.maxInexactHomLength)
            val maxInexactHomLengthShortDel = defaultIntValue(cmd, MAX_INEXACT_HOM_LENGTH_SHORT_DEL_OPTION, defaultConfig.maxInexactHomLengthShortDel)
            val minLength = defaultIntValue(cmd, MIN_LENGTH_OPTION, defaultConfig.minLength)

            return GripssFilterConfig(
                    maxNormalSupportProportion,
                    minNormalCoverage,
                    minTumorAF,
                    maxShortStrandBias,
                    minQualBreakEnd,
                    minQualBreakPoint,
                    maxHomLengthShortInversion,
                    maxInexactHomLength,
                    maxInexactHomLengthShortDel,
                    minLength)
        }

        private fun default(): GripssFilterConfig {
            return GripssFilterConfig(
                    0.03,
                    8,
                    0.005,
                    0.95,
                    1000,
                    350,
                    6,
                    50,
                    5,
                    32)
        }
    }
}