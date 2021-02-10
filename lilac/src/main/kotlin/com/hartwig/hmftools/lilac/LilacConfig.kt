package com.hartwig.hmftools.lilac

import com.hartwig.hmftools.common.cli.Configs
import com.hartwig.hmftools.lilac.hla.HlaAllele
import org.apache.commons.cli.CommandLine
import org.apache.commons.cli.Option
import org.apache.commons.cli.Options
import org.apache.commons.cli.ParseException
import java.io.File
import java.io.IOException

const val SAMPLE = "sample"
const val INPUT_BAM_OPTION = "sample_bam"
const val RESOURCE_DIR_OPTION = "resource_dir"
const val OUTPUT_DIR_OPTION = "output_dir"
const val REF_GENOME_OPTION = "ref_genome"
const val MIN_BASE_QUAL = "min_base_qual"
const val MIN_EVIDENCE = "min_evidence"
const val THREADS = "threads"
const val MIN_FRAGMENTS_PER_ALLELE = "min_fragments_per_allele"
const val MIN_FRAGMENTS_TO_REMOVE_SINGLE = "min_fragments_to_remove_single"
const val MIN_CONFIRMED_UNIQUE_COVERAGE = "min_confirmed_unique_coverage"
const val EXPECTED_ALLELES = "expected_alleles"

data class LilacConfig(
        val sample: String,
        val inputBam: String,
        val resourceDir: String,
        val outputDir: String,
        val refGenome: String,
        val minBaseQual: Int,
        val minEvidence: Int,
        val minFragmentsPerAllele: Int,
        val minFragmentsToRemoveSingle: Int,
        val minConfirmedUniqueCoverage: Int,
        val threads: Int,
        val expectedAlleles: List<HlaAllele>
) {

    val outputFilePrefix = "${outputDir}/$sample"

    companion object {

        @Throws(ParseException::class, IOException::class)
        fun createConfig(cmd: CommandLine): LilacConfig {
            val sample = cmd.getOptionValue(SAMPLE)
            val inputBam = cmd.requiredFile(INPUT_BAM_OPTION)
            val resourceDir = cmd.requiredDir(RESOURCE_DIR_OPTION)
            val outputDir = cmd.requiredDir(OUTPUT_DIR_OPTION)
            val defaultConfig = default()

            val refGenome = cmd.optionalFile(REF_GENOME_OPTION, "")
            val minBaseQual = Configs.defaultIntValue(cmd, MIN_BASE_QUAL, defaultConfig.minBaseQual)
            val minEvidence = Configs.defaultIntValue(cmd, MIN_EVIDENCE, defaultConfig.minEvidence)
            val minFragmentsPerAllele = Configs.defaultIntValue(cmd, MIN_FRAGMENTS_PER_ALLELE, defaultConfig.minFragmentsPerAllele)
            val minFragmentsToRemoveSingle = Configs.defaultIntValue(cmd, MIN_FRAGMENTS_TO_REMOVE_SINGLE, defaultConfig.minFragmentsToRemoveSingle)
            val minConfirmedUniqueCoverage = Configs.defaultIntValue(cmd, MIN_CONFIRMED_UNIQUE_COVERAGE, defaultConfig.minConfirmedUniqueCoverage)
            val threads = Configs.defaultIntValue(cmd, THREADS, defaultConfig.threads)
            val expectedAlleleString = Configs.defaultStringValue(cmd, EXPECTED_ALLELES, "")
            val expectedAlleles = if (expectedAlleleString.isNotEmpty()) {
                expectedAlleleString.split(",").map { HlaAllele(it) }
            } else {
                listOf()
            }


            return LilacConfig(
                    sample,
                    inputBam,
                    resourceDir,
                    outputDir,
                    refGenome,
                    minBaseQual,
                    minEvidence,
                    minFragmentsPerAllele,
                    minFragmentsToRemoveSingle,
                    minConfirmedUniqueCoverage,
                    threads,
                    expectedAlleles)
        }

        private fun default(): LilacConfig {
            return LilacConfig(
                    "",
                    "",
                    "",
                    "",
                    "",
                    30,
                    3,
                    6,
                    40,
                    10,
                    1,
                    listOf())
        }

        fun createOptions(): Options {
            val options = Options()
            options.addOption(requiredOption(SAMPLE, "Name of sample"))
            options.addOption(requiredOption(INPUT_BAM_OPTION, "Path to input bam"))
            options.addOption(requiredOption(RESOURCE_DIR_OPTION, "Path to resource files"))
            options.addOption(requiredOption(OUTPUT_DIR_OPTION, "Path to output"))
            options.addOption(optional(REF_GENOME_OPTION, "REF_GENOME_OPTION"))
            options.addOption(optional(MIN_BASE_QUAL, "MIN_BASE_QUAL"))
            options.addOption(optional(MIN_EVIDENCE, "MIN_EVIDENCE"))
            options.addOption(optional(MIN_FRAGMENTS_PER_ALLELE, "MIN_FRAGMENTS_PER_ALLELE"))
            options.addOption(optional(MIN_FRAGMENTS_TO_REMOVE_SINGLE, "MIN_FRAGMENTS_TO_REMOVE_SINGLE"))
            options.addOption(optional(MIN_CONFIRMED_UNIQUE_COVERAGE, "MIN_CONFIRMED_UNIQUE_COVERAGE"))
            options.addOption(optional(THREADS, "Number of threads"))
            options.addOption(optional(EXPECTED_ALLELES, "Common separated expected alleles"))
            return options
        }

        @Throws(IOException::class)
        internal fun CommandLine.requiredFile(argument: String): String {
            val result = this.getOptionValue(argument)
            if (!File(result).exists()) {
                throw IOException("Unable to read file $result")
            }

            return result
        }

        @Throws(IOException::class)
        internal fun CommandLine.optionalFile(argument: String, default: String): String {
            if (this.hasOption(argument)) {
                return this.requiredFile(argument)
            }

            return default
        }

        @Throws(IOException::class)
        internal fun CommandLine.requiredDir(argument: String): String {
            val result = this.getOptionValue(argument)
            val dir = File(result)
            if (!dir.exists() && !dir.mkdirs()) {
                throw IOException("Unable to create director $result")
            }

            return result
        }

        private fun optional(argument: String, description: String): Option {
            val result = Option(argument, true, description)
            result.isRequired = false
            return result
        }


        private fun requiredOption(argument: String, description: String): Option {
            val result = Option(argument, true, description)
            result.isRequired = true
            return result
        }

    }

}