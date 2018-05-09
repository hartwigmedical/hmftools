package com.hartwig.hmftools.actionabilityAnalyzer

import com.hartwig.hmftools.patientdb.dao.DatabaseAccess
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import org.apache.logging.log4j.LogManager
import java.io.FileWriter
import kotlin.streams.asSequence

private val logger = LogManager.getLogger("ActionabilityAnalyzerApplication")
private val user = ""
private val password = ""
private val databaseUrl = "jdbc:mysql://"

fun main(args: Array<String>) {
    val actionabilityAnalyzer = ActionabilityAnalyzer("actionableVariantsFinal")
    val dbAccess = DatabaseAccess(user, password, databaseUrl)
    val printer = createPrinter()
    logger.info("Start")
    queryDatabase(dbAccess, printer, actionabilityAnalyzer)
    logger.info("Done.")
}

private fun queryDatabase(dbAccess: DatabaseAccess, printer: CSVPrinter, actionabilityAnalyzer: ActionabilityAnalyzer) {
    dbAccess.potentiallyActionableVariants().use {
        it.asSequence().forEachIndexed { index, variant ->
            logProgress(index)
            actionabilityAnalyzer.actionabilityForVariant(variant).forEach { printer.printRecord(it.record) }
            printer.flush()
        }
    }
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
