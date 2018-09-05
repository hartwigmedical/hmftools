package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.extensions.cli.createCommandLine
import com.hartwig.hmftools.extensions.cli.options.HmfOptions
import com.hartwig.hmftools.extensions.cli.options.strings.RequiredInputOption
import com.hartwig.hmftools.extensions.csv.CsvWriter
import com.hartwig.hmftools.knowledgebaseimporter.DB_PASSWORD
import com.hartwig.hmftools.knowledgebaseimporter.DB_USER
import com.hartwig.hmftools.knowledgebaseimporter.HMFPATIENTS_DB
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess
import com.hartwig.hmftools.patientdb.data.PotentialActionableCNV
import com.hartwig.hmftools.patientdb.data.PotentialActionableFusion
import org.apache.logging.log4j.LogManager
import java.io.File
import java.util.stream.Stream
import kotlin.streams.asSequence
import kotlin.streams.toList

private val logger = LogManager.getLogger("ActionabilityAnalyzerApplication")
const val SAMPLE_ID = "sample_ID"
const val OUTPUT_DIRECTORY = "output_dir"
const val KNOWLEDGEBASE = "knowledgebase"

fun main(args: Array<String>) {
    logger.info("Start processing actionability")
    val cmd = createOptions().createCommandLine("actionability-analyzer", args)
    val sampleId = cmd.getOptionValue(SAMPLE_ID)
    val outputDir = cmd.getOptionValue(OUTPUT_DIRECTORY)
    val knowledgebase = cmd.getOptionValue(KNOWLEDGEBASE)
    val user = cmd.getOptionValue(DB_USER)
    val password = cmd.getOptionValue(DB_PASSWORD)
    val databaseUrl = "jdbc:${cmd.getOptionValue(HMFPATIENTS_DB)}"

    val actionableVariants = "$knowledgebase/actionableVariants.tsv"
    val actionableFusionPairs = "$knowledgebase/actionableFusionPairs.tsv"
    val actionablePromiscuousFive = "$knowledgebase/actionablePromiscuousFive.tsv"
    val actionablePromiscuousThree = "$knowledgebase/actionablePromiscuousThree.tsv"
    val actionableCNVs = "$knowledgebase/actionableCNVs.tsv"
    val actionableGenomicRanges = "$knowledgebase/actionableRanges.tsv"
    val cancerTypes = "$knowledgebase/knowledgebaseCancerTypes.tsv"

    val dbAccess = DatabaseAccess(user, password, databaseUrl)
    val sampleToAnalyze = readSample(sampleId, dbAccess)
    val actionabilityAnalyzer = ActionabilityAnalyzer(sampleToAnalyze, actionableVariants, actionableFusionPairs,
                                                        actionablePromiscuousFive, actionablePromiscuousThree, actionableCNVs, cancerTypes,
                                                        actionableGenomicRanges)
    queryDatabase(outputDir, dbAccess, sampleToAnalyze, actionabilityAnalyzer)
    logger.info("Done processing actionability")
}

private fun createOptions(): HmfOptions {
    val options = HmfOptions()
    options.add(RequiredInputOption(SAMPLE_ID, "sample_ID"))
    options.add(RequiredInputOption(OUTPUT_DIRECTORY, "output_dir"))
    options.add(RequiredInputOption(KNOWLEDGEBASE, "knowledgebase"))
    options.add(RequiredInputOption(HMFPATIENTS_DB, "hmfpatients db url"))
    options.add(RequiredInputOption(DB_USER, "db user"))
    options.add(RequiredInputOption(DB_PASSWORD, "db password"))
    return options
}

private fun readSample(sampleId: String, dbAccess: DatabaseAccess) : Map<String, String> =
        if (sampleId.isNotEmpty() && sampleId.isNotBlank()) dbAccess.sampleAndTumorLocation(sampleId).toList().associate { Pair(it.key, it.value) }
        else mapOf(sampleId to "It is no tumor sample")

private fun queryDatabase(outputDir: String, dbAccess: DatabaseAccess, sampleToAnalyze: Map<String, String>, actionabilityAnalyzer: ActionabilityAnalyzer) {
    val records = mutableListOf<ActionabilityOutput>()
    potentiallyActionableCNVs(dbAccess, sampleToAnalyze).use {
        val cnvRecords = it.asSequence().flatMap { actionabilityAnalyzer.actionabilityForCNV(it).asSequence() }.toList()
        records.addAll(cnvRecords)
    }
    potentiallyActionableFusions(dbAccess, sampleToAnalyze).use {
        val fusionRecords = it.asSequence().flatMap { actionabilityAnalyzer.actionabilityForFusion(it).asSequence() }.toList()
        records.addAll(fusionRecords)
    }
    logger.info("Done. Writing results to file...")
    CsvWriter.writeTSV(records, "$outputDir${File.separator}actionableVariantsPerSample.tsv")
}

private fun potentiallyActionableCNVs(dbAccess: DatabaseAccess, sampleToAnalyze: Map<String, String>): Stream<PotentialActionableCNV> {
    logger.info("Querying actionable cnvs.")
    return if (false) dbAccess.allPotentiallyActionableCNVs()
    else dbAccess.potentiallyActionableCNVs(sampleToAnalyze.keys)
}

private fun potentiallyActionableFusions(dbAccess: DatabaseAccess,
                                         sampleToAnalyze: Map<String, String>): Stream<PotentialActionableFusion> {
    logger.info("Querying actionable fusions.")
    return if (false) dbAccess.allPotentiallyActionableFusions()
    else dbAccess.potentiallyActionableFusions(sampleToAnalyze.keys)
}