package com.hartwig.hmftools.lilac

import com.hartwig.hmftools.common.cli.Configs
import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess
import org.apache.commons.cli.CommandLine
import org.apache.commons.cli.Option
import org.apache.commons.cli.Options
import org.apache.commons.cli.ParseException
import org.apache.logging.log4j.LogManager
import java.io.File
import java.io.IOException

const val SAMPLE = "sample"
const val REFERENCE_BAM_OPTION = "reference_bam"
const val TUMOR_BAM_OPTION = "tumor_bam"
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
const val GENE_COPY_NUMBER = "gene_copy_number"
const val SOMATIC_VCF = "somatic_vcf"
const val MAX_DISTANCE_FROM_TOP_SCORE = "max_distance_from_top_score"
const val DEBUG_PHASING = "debug_phasing"

data class LilacConfig(
        val sample: String,
        val referenceBam: String,
        val tumorBam: String,
        val resourceDir: String,
        val outputDir: String,
        val refGenome: String,
        val minBaseQual: Int,
        val minEvidence: Int,
        val minFragmentsPerAllele: Int,
        val minFragmentsToRemoveSingle: Int,
        val minConfirmedUniqueCoverage: Int,
        val threads: Int,
        val maxDistanceFromTopScore: Int,
        val geneCopyNumberFile: String,
        val somaticVcf: String,
        val debugPhasing: Boolean,
        val expectedAlleles: List<HlaAllele>,
        val commonAlleles: List<HlaAllele>
) {

    val outputFilePrefix = "${outputDir}/$sample"

    companion object {

        val logger = LogManager.getLogger(this::class.java)

        @Throws(ParseException::class, IOException::class)
        fun createConfig(cmd: CommandLine): LilacConfig {
            val sample = cmd.getOptionValue(SAMPLE)
            val referenceBam = cmd.requiredFile(REFERENCE_BAM_OPTION)
            val tumorBam = cmd.optionalFile(TUMOR_BAM_OPTION, "")
            val geneCopyNumberFile = cmd.optionalFile(GENE_COPY_NUMBER, "")
            val somaticVcf = cmd.optionalFile(SOMATIC_VCF, "")
            val resourceDir = cmd.requiredDir(RESOURCE_DIR_OPTION)
            val outputDir = cmd.requiredDir(OUTPUT_DIR_OPTION)
            val defaultConfig = default()

            val debugPhasing = cmd.hasOption(DEBUG_PHASING)
            val refGenome = cmd.optionalFile(REF_GENOME_OPTION, "")
            val minBaseQual = Configs.defaultIntValue(cmd, MIN_BASE_QUAL, defaultConfig.minBaseQual)
            val minEvidence = Configs.defaultIntValue(cmd, MIN_EVIDENCE, defaultConfig.minEvidence)
            val minFragmentsPerAllele = Configs.defaultIntValue(cmd, MIN_FRAGMENTS_PER_ALLELE, defaultConfig.minFragmentsPerAllele)
            val minFragmentsToRemoveSingle = Configs.defaultIntValue(cmd, MIN_FRAGMENTS_TO_REMOVE_SINGLE, defaultConfig.minFragmentsToRemoveSingle)
            val minConfirmedUniqueCoverage = Configs.defaultIntValue(cmd, MIN_CONFIRMED_UNIQUE_COVERAGE, defaultConfig.minConfirmedUniqueCoverage)
            val threads = Configs.defaultIntValue(cmd, THREADS, defaultConfig.threads)
            val maxDistanceFromTopScore = Configs.defaultIntValue(cmd, MAX_DISTANCE_FROM_TOP_SCORE, defaultConfig.maxDistanceFromTopScore)
            val expectedAlleles = cmd.expectedAlleles(EXPECTED_ALLELES);
            val commonAlleles = LilacConfig::class.java.getResource("/alleles/common.txt")
                    .readText()
                    .split("\n")
                    .map { HlaAllele(it) }
                    .map { it.asFourDigit() }

            return LilacConfig(
                    sample,
                    referenceBam,
                    tumorBam,
                    resourceDir,
                    outputDir,
                    refGenome,
                    minBaseQual,
                    minEvidence,
                    minFragmentsPerAllele,
                    minFragmentsToRemoveSingle,
                    minConfirmedUniqueCoverage,
                    threads,
                    maxDistanceFromTopScore,
                    geneCopyNumberFile,
                    somaticVcf,
                    debugPhasing,
                    expectedAlleles,
                    commonAlleles)
        }

        private fun default(): LilacConfig {
            return LilacConfig(
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    30,
                    2,
                    6,
                    40,
                    10,
                    1,
                    3,
                    "",
                    "",
                    false,
                    listOf(),
                    listOf())
        }

        fun createOptions(): Options {
            val options = Options()
            options.addOption(requiredOption(SAMPLE, "Name of sample"))
            options.addOption(requiredOption(REFERENCE_BAM_OPTION, "Path to reference/normal bam"))
            options.addOption(optional(TUMOR_BAM_OPTION, "Path to tumor bam"))
            options.addOption(requiredOption(RESOURCE_DIR_OPTION, "Path to resource files"))
            options.addOption(requiredOption(OUTPUT_DIR_OPTION, "Path to output"))
            options.addOption(optional(REF_GENOME_OPTION, "Optional path to reference genome fasta file"))
            options.addOption(optional(MIN_BASE_QUAL, "MIN_BASE_QUAL"))
            options.addOption(optional(MIN_EVIDENCE, "MIN_EVIDENCE"))
            options.addOption(optional(MIN_FRAGMENTS_PER_ALLELE, "MIN_FRAGMENTS_PER_ALLELE"))
            options.addOption(optional(MIN_FRAGMENTS_TO_REMOVE_SINGLE, "MIN_FRAGMENTS_TO_REMOVE_SINGLE"))
            options.addOption(optional(MIN_CONFIRMED_UNIQUE_COVERAGE, "MIN_CONFIRMED_UNIQUE_COVERAGE"))
            options.addOption(optional(MAX_DISTANCE_FROM_TOP_SCORE, "Max distance from top score"))
            options.addOption(optional(THREADS, "Number of threads"))
            options.addOption(optional(EXPECTED_ALLELES, "Comma separated expected alleles"))
            options.addOption(optional(GENE_COPY_NUMBER, "Path to gene copy number file"))
            options.addOption(optional(SOMATIC_VCF, "Path to somatic VCF"))
            options.addOption(Option(DEBUG_PHASING, false, "More detailed logging of phasing"))
            DatabaseAccess.addDatabaseCmdLineArgs(options)
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

        private fun CommandLine.expectedAlleles(opt: String): List<HlaAllele> {
            if (this.hasOption(opt)) {
                val result = this.getOptionValue(opt)!!.split(",").map { HlaAllele(it).asFourDigit() }
                logger.info("Using non default value {} for parameter {}", result.joinToString(","), opt)
                return result
            }
            return listOf()
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