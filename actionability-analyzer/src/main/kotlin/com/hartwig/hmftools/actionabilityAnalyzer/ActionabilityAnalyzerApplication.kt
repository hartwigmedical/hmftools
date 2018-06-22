package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.extensions.csv.CsvWriter
import com.hartwig.hmftools.knowledgebaseimporter.readCSVRecords
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess
import com.hartwig.hmftools.patientdb.data.PotentialActionableCNV
import com.hartwig.hmftools.patientdb.data.PotentialActionableFusion
import com.hartwig.hmftools.patientdb.data.PotentialActionableVariant
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
    logger.info("Querying actionable variants.")
    val records = mutableListOf<ActionabilityOutput>()
    potentiallyActionableVariants(dbAccess, samplesToAnalyze).use {
        val variantRecords = it.asSequence().mapIndexed { index, variant ->
            logProgress(index)
            actionabilityAnalyzer.actionabilityForVariant(variant) + actionabilityAnalyzer.rangeActionabilityForVariant(variant)
        }.flatten().toList()
        records.addAll(variantRecords)
    }
    logger.info("Querying actionable cnvs.")
    potentiallyActionableCNVs(dbAccess, samplesToAnalyze).use {
        val cnvRecords = it.asSequence().flatMap { actionabilityAnalyzer.actionabilityForCNV(it).asSequence() }.toList()
        records.addAll(cnvRecords)
    }
    logger.info("Querying actionable fusions.")
    potentiallyActionableFusions(dbAccess, samplesToAnalyze).use {
        val fusionRecords = it.asSequence().flatMap { actionabilityAnalyzer.actionabilityForFusion(it).asSequence() }.toList()
        records.addAll(fusionRecords)
    }
    logger.info("Done. Writing results to file...")
    CsvWriter.writeTSV(records, "actionableVariantsPerSample.tsv")
    logger.info("Done.")
}

private fun logProgress(index: Int) {
    if (index % 1000000 == 0) {
        logger.info("Processed $index records...")
    }
}

private fun readSamples(dbAccess: DatabaseAccess): Map<String, String> =
        if (allDbSamples) dbAccess.allSamplesAndTumorLocations().toList().associate { Pair(it.key, it.value) }
        else readCSVRecords(cohortCsvLocation) { Pair(it[0], it[1]) }.toMap()

private fun potentiallyActionableVariants(dbAccess: DatabaseAccess,
                                          samplesToAnalyze: Map<String, String>): Stream<PotentialActionableVariant> =
        if (allDbSamples) dbAccess.allPotentiallyActionableVariants()
        else dbAccess.potentiallyActionableVariants(samplesToAnalyze.keys)

private fun potentiallyActionableCNVs(dbAccess: DatabaseAccess,
                                      samplesToAnalyze: Map<String, String>): Stream<PotentialActionableCNV> =
        if (allDbSamples) dbAccess.allPotentiallyActionableCNVs()
        else dbAccess.potentiallyActionableCNVs(samplesToAnalyze.keys)

private fun potentiallyActionableFusions(dbAccess: DatabaseAccess,
                                         samplesToAnalyze: Map<String, String>): Stream<PotentialActionableFusion> =
        if (allDbSamples) dbAccess.allPotentiallyActionableFusions()
        else dbAccess.potentiallyActionableFusions(samplesToAnalyze.keys)