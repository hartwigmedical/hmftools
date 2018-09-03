package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.extensions.cli.createCommandLine
import com.hartwig.hmftools.extensions.cli.options.HmfOptions
import com.hartwig.hmftools.extensions.cli.options.strings.RequiredInputOption
import com.hartwig.hmftools.extensions.csv.CsvReader
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
private val actionableVariants = "actionableVariants"
private val actionableFusionPairs = "actionableFusionPairs"
private val actionablePromiscuousFive = "actionablePromiscuousFive"
private val actionablePromiscuousThree = "actionablePromiscuousThree"
private val actionableCNVs = "actionableCNVs"
private val actionableGenomicRanges = "actionableRanges"
private val cancerTypes = "knowledgebaseCancerTypes"
private val allDbSamples = false
private val cohortCsvLocation = "cohort.csv"
private val cohortMutations = "cohortMutations.tsv"
const val SAMPLE_ID = "sample_ID"
const val OUTPUT_DIRECTORY = "/data/common/dbs/knowledgebases/actionability"

fun main(args: Array<String>) {
    logger.info("Start")
    val cmd = createOptions().createCommandLine("actionability-analyzer", args)
    val sampleId = cmd.getOptionValue(SAMPLE_ID)
    logger.info(sampleId)
    val user = cmd.getOptionValue(DB_USER)
    val password = cmd.getOptionValue(DB_PASSWORD)
    val databaseUrl = "jdbc:${cmd.getOptionValue(HMFPATIENTS_DB)}"

    val dbAccess = DatabaseAccess(user, password, databaseUrl)
    val samplesToAnalyze = readSamples(sampleId, dbAccess)
    logger.info("samplesToAnalyze" + samplesToAnalyze)
    //val actionabilityAnalyzer = ActionabilityAnalyzer(samplesToAnalyze.toString(), actionableVariants, actionableFusionPairs,
      //                                                  actionablePromiscuousFive, actionablePromiscuousThree, actionableCNVs, cancerTypes,
        //                                                actionableGenomicRanges)
    //logger.info("actionabilityAnalyzer" + actionabilityAnalyzer)
    //   queryDatabase(outputDir, dbAccess, samplesToAnalyze, actionabilityAnalyzer)
    logger.info("Done.")
}

private fun createOptions(): HmfOptions {
    val options = HmfOptions()
    options.add(RequiredInputOption(SAMPLE_ID, "sample_ID"))
    options.add(RequiredInputOption(HMFPATIENTS_DB, "hmfpatients db url"))
    options.add(RequiredInputOption(DB_USER, "db user"))
    options.add(RequiredInputOption(DB_PASSWORD, "db password"))
    return options
}

private fun readSamples(sampleId: String, dbAccess: DatabaseAccess) =
        if (sampleId.isNotEmpty() && sampleId.isNotBlank()) dbAccess.allSamplesAndTumorLocations(sampleId).toList().associate { Pair(it.key, it.value) }
        else logger.info("no sampleId")

private fun queryDatabase(outputDir: String, dbAccess: DatabaseAccess, samplesToAnalyze: Map<String, String>, actionabilityAnalyzer: ActionabilityAnalyzer) {
    val records = mutableListOf<ActionabilityOutput>()
    cohortMutations(samplesToAnalyze).forEach { mutation ->
        val variantRecords = actionabilityAnalyzer.actionabilityForVariant(mutation) +
                actionabilityAnalyzer.rangeActionabilityForVariant(mutation)
        records.addAll(variantRecords)
    }
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
    logger.info("Done.")
}

private fun cohortMutations(samplesToAnalyze: Map<String, String>): List<CohortMutation> {
    logger.info("Looking up cohort mutations")
    return CsvReader.readTSV<CohortMutation>(cohortMutations).filter { samplesToAnalyze.containsKey(it.sampleId) }
            .flatMap { mutation -> mutation.alt.split(",").map { mutation.copy(alt = it) } }
}

private fun potentiallyActionableCNVs(dbAccess: DatabaseAccess, samplesToAnalyze: Map<String, String>): Stream<PotentialActionableCNV> {
    logger.info("Querying actionable cnvs.")
    return if (allDbSamples) dbAccess.allPotentiallyActionableCNVs()
    else dbAccess.potentiallyActionableCNVs(samplesToAnalyze.keys)
}

private fun potentiallyActionableFusions(dbAccess: DatabaseAccess,
                                         samplesToAnalyze: Map<String, String>): Stream<PotentialActionableFusion> {
    logger.info("Querying actionable fusions.")
    return if (allDbSamples) dbAccess.allPotentiallyActionableFusions()
    else dbAccess.potentiallyActionableFusions(samplesToAnalyze.keys)
}