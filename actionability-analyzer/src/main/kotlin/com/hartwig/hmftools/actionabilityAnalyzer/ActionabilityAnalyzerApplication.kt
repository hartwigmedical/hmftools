package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.knowledgebaseimporter.readCSVRecords
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import org.apache.logging.log4j.LogManager
import java.io.FileWriter
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
private val cancerTypes = "knowledgebaseCancerTypes"
private val allDbSamples = false
private val cohortCsvLocation = "cohort.csv"

fun main(args: Array<String>) {
    val dbAccess = DatabaseAccess(user, password, databaseUrl)
    val samplesToAnalyze = readSamples(dbAccess)
    val actionabilityAnalyzer = ActionabilityAnalyzer(samplesToAnalyze, actionableVariants, actionableFusionPairs,
                                                      actionablePromiscuousFive, actionablePromisucousThree, actionableCNVs, cancerTypes)
    val printer = createPrinter()
    logger.info("Start")
    queryDatabase(dbAccess, printer, actionabilityAnalyzer)
    logger.info("Done.")
}

private fun queryDatabase(dbAccess: DatabaseAccess, printer: CSVPrinter, actionabilityAnalyzer: ActionabilityAnalyzer) {
    logger.info("Querying actionable variants.")
    dbAccess.potentiallyActionableVariants().use {
        it.asSequence().forEachIndexed { index, variant ->
            logProgress(index)
            val records = actionabilityAnalyzer.actionabilityForVariant(variant).map { it.record }
            printer.printRecords(records)
            printer.flush()
        }
    }
    logger.info("Done writing actionable variants.")
    logger.info("Querying actionable cnvs.")
    dbAccess.potentiallyActionableCNVs().use {
        it.asSequence().forEach { cnv ->
            val records = actionabilityAnalyzer.actionabilityForCNV(cnv).map { it.record }
            printer.printRecords(records)
            printer.flush()
        }
    }
    logger.info("Done writing actionable cnvs.")
    logger.info("Querying actionable fusions.")
    dbAccess.potentiallyActionableFusions().use {
        it.asSequence().forEach { fusion ->
            val records = actionabilityAnalyzer.actionabilityForFusion(fusion).map { it.record }
            printer.printRecords(records)
            printer.flush()
        }
    }
    logger.info("Done writing actionable fusions.")
    printer.close()
}

private fun logProgress(index: Int) {
    if (index % 1000000 == 0) {
        logger.info("Processed $index records...")
    }
}

private fun createPrinter(): CSVPrinter {
    val format = CSVFormat.TDF.withHeader(*ActionabilityOutput.header.toTypedArray()).withNullString("NULL")
    return CSVPrinter(FileWriter("actionableVariantsPerSample.tsv"), format)
}

private fun readSamples(dbAccess: DatabaseAccess): Map<String, String> =
        if (allDbSamples) dbAccess.allSamplesAndTumorLocations().toList().associate { Pair(it.key, it.value) }
        else readCSVRecords(cohortCsvLocation) { Pair(it[0], it[1]) }.toMap()
