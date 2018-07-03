package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.extensions.csv.CsvReader
import com.hartwig.hmftools.extensions.csv.CsvWriter
import com.hartwig.hmftools.knowledgebaseimporter.readCSVRecords
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess
import com.hartwig.hmftools.patientdb.data.PotentialActionableCNV
import com.hartwig.hmftools.patientdb.data.PotentialActionableFusion
import org.apache.logging.log4j.LogManager
import java.util.stream.Stream
import kotlin.streams.asSequence
import kotlin.streams.toList

private val logger = LogManager.getLogger("ActionabilityAnalyzerApplication")
private val user = ""
private val password = ""
private val databaseUrl = "jdbc:mysql://"
private val actionableVariants = "actionableVariants"
private val actionableFusionPairs = "actionableFusionPairs"
private val actionablePromiscuousFive = "actionablePromiscuousFive"
private val actionablePromisucousThree = "actionablePromiscuousThree"
private val actionableCNVs = "actionableCNVs"
private val actionableGenomicRanges = "actionableRanges"
private val cancerTypes = "knowledgebaseCancerTypes"
private val allDbSamples = false
private val cohortCsvLocation = "cohort.csv"
private val cohortMutations = "cohortMutations.tsv"

fun main(args: Array<String>) {
    val dbAccess = DatabaseAccess(user, password, databaseUrl)
    val samplesToAnalyze = readSamples(dbAccess)
    val actionabilityAnalyzer = ActionabilityAnalyzer(samplesToAnalyze, actionableVariants, actionableFusionPairs,
                                                      actionablePromiscuousFive, actionablePromisucousThree, actionableCNVs, cancerTypes,
                                                      actionableGenomicRanges)
    logger.info("Start")
    queryDatabase(dbAccess, samplesToAnalyze, actionabilityAnalyzer)
    logger.info("Done.")
}

private fun queryDatabase(dbAccess: DatabaseAccess, samplesToAnalyze: Map<String, String>, actionabilityAnalyzer: ActionabilityAnalyzer) {
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
    CsvWriter.writeTSV(records, "actionableVariantsPerSample.tsv")
    logger.info("Done.")
}

private fun readSamples(dbAccess: DatabaseAccess): Map<String, String> =
        if (allDbSamples) dbAccess.allSamplesAndTumorLocations().toList().associate { Pair(it.key, it.value) }
        else readCSVRecords(cohortCsvLocation) { Pair(it[0], it[1]) }.toMap()

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
