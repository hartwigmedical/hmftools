package com.hartwig.hmftools.actionabilityAnalyzer

import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVPrinter
import org.apache.logging.log4j.LogManager
import java.io.FileWriter

private val logger = LogManager.getLogger("ActionabilityAnalyzerApplication")

fun main(args: Array<String>) {
    
}

private fun createPrinter(): CSVPrinter {
    val format = CSVFormat.TDF.withHeader(*ActionabilityOutput.header.toTypedArray()).withNullString("NULL")
    return CSVPrinter(FileWriter("actionableVariantsPerSample.tsv"), format)
}
