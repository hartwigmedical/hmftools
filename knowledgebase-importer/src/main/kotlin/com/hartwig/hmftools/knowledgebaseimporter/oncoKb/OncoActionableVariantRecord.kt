package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import org.apache.commons.csv.CSVRecord

data class OncoActionableVariantRecord(private val csvRecord: CSVRecord) {
    val transcript: String = csvRecord["Isoform"]
    val gene: String = csvRecord["Gene"]
    val alteration: String = csvRecord["Alteration"]
    val cancerType: String = csvRecord["Cancer Type"]
    val drugs: List<String> = csvRecord["Drugs(s)"].split(",").map { it.trim() }
    val level: String = if (csvRecord["Level"].startsWith("R")) {
        csvRecord["Level"].substring(1)
    } else {
        csvRecord["Level"]
    }
    val significance = if (csvRecord["Level"].startsWith("R")) {
        "resistance"
    } else {
        "sensitivity"
    }
}
