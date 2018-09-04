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
private val actionableVariants = "/data/common/dbs/knowledgebases/output_180822/actionableVariants.tsv"
private val actionableFusionPairs = "/data/common/dbs/knowledgebases/output_180822/actionableFusionPairs.tsv"
private val actionablePromiscuousFive = "/data/common/dbs/knowledgebases/output_180822/actionablePromiscuousFive.tsv"
private val actionablePromiscuousThree = "/data/common/dbs/knowledgebases/output_180822/actionablePromiscuousThree.tsv"
private val actionableCNVs = "/data/common/dbs/knowledgebases/output_180822/actionableCNVs.tsv"
private val actionableGenomicRanges = "/data/common/dbs/knowledgebases/output_180822/actionableRanges.tsv"
private val cancerTypes = "/data/common/dbs/knowledgebases/output_180822/knowledgebaseCancerTypes.tsv"
const val SAMPLE_ID = "sample_ID"
const val OUTPUT_DIRECTORY = "/data/common/dbs/knowledgebases/actionability/"

fun main(args: Array<String>) {
    logger.info("Start processing actionability")
    val cmd = createOptions().createCommandLine("actionability-analyzer", args)
    val sampleId = cmd.getOptionValue(SAMPLE_ID)
    logger.info(sampleId)
    val user = cmd.getOptionValue(DB_USER)
    val password = cmd.getOptionValue(DB_PASSWORD)
    val databaseUrl = "jdbc:${cmd.getOptionValue(HMFPATIENTS_DB)}"

    val dbAccess = DatabaseAccess(user, password, databaseUrl)
    val samplesToAnalyze = readSamples(sampleId, dbAccess)
    val actionabilityAnalyzer = ActionabilityAnalyzer(samplesToAnalyze, actionableVariants, actionableFusionPairs,
                                                        actionablePromiscuousFive, actionablePromiscuousThree, actionableCNVs, cancerTypes,
                                                        actionableGenomicRanges)
    queryDatabase(OUTPUT_DIRECTORY, dbAccess, samplesToAnalyze, actionabilityAnalyzer)
    logger.info("Done processing actionability")
}

private fun createOptions(): HmfOptions {
    val options = HmfOptions()
    options.add(RequiredInputOption(SAMPLE_ID, "sample_ID"))
    options.add(RequiredInputOption(HMFPATIENTS_DB, "hmfpatients db url"))
    options.add(RequiredInputOption(DB_USER, "db user"))
    options.add(RequiredInputOption(DB_PASSWORD, "db password"))
    return options
}

private fun readSamples(sampleId: String, dbAccess: DatabaseAccess) : Map<String, String> =
        if (sampleId.isNotEmpty() && sampleId.isNotBlank()) dbAccess.allSamplesAndTumorLocations(sampleId).toList().associate { Pair(it.key, it.value) }
        else mapOf(sampleId to "It is no tumor sample")

private fun queryDatabase(outputDir: String, dbAccess: DatabaseAccess, samplesToAnalyze: Map<String, String>, actionabilityAnalyzer: ActionabilityAnalyzer) {
    val records = mutableListOf<ActionabilityOutput>()
    potentiallyActionableCNVs(dbAccess, samplesToAnalyze).use {
        val cnvRecords = it.asSequence().flatMap { actionabilityAnalyzer.actionabilityForCNV(it).asSequence() }.toList()
        records.addAll(cnvRecords)
    }
    potentiallyActionableFusions(dbAccess, samplesToAnalyze).use {
        val fusionRecords = it.asSequence().flatMap { actionabilityAnalyzer.actionabilityForFusion(it).asSequence() }.toList()
        records.addAll(fusionRecords)
    }
    logger.info("Done. Writing results to file...")
    CsvWriter.writeTSV(records, "$outputDir${File.separator}actionableVariantsPerSample.tsv")
}

private fun potentiallyActionableCNVs(dbAccess: DatabaseAccess, samplesToAnalyze: Map<String, String>): Stream<PotentialActionableCNV> {
    logger.info("Querying actionable cnvs.")
    return if (false) dbAccess.allPotentiallyActionableCNVs()
    else dbAccess.potentiallyActionableCNVs(samplesToAnalyze.keys)
}

private fun potentiallyActionableFusions(dbAccess: DatabaseAccess,
                                         samplesToAnalyze: Map<String, String>): Stream<PotentialActionableFusion> {
    logger.info("Querying actionable fusions.")
    return if (false) dbAccess.allPotentiallyActionableFusions()
    else dbAccess.potentiallyActionableFusions(samplesToAnalyze.keys)
}