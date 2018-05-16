package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfLevel
import org.apache.commons.csv.CSVRecord

data class OncoActionableVariantRecord(private val csvRecord: CSVRecord) {
    val transcript: String = csvRecord["Isoform"]
    val gene: String = csvRecord["Gene"]
    val alteration: String = csvRecord["Alteration"]
    val cancerType: String = csvRecord["Cancer Type"]
    val drugs: List<String> = csvRecord["Drugs(s)"].split(",").map { it.trim() }
    val level: String = readLevel(csvRecord["Level"])
    val significance = if (csvRecord["Level"].startsWith("R")) "resistance" else "sensitivity"
    val actionabilityItems: List<Actionability> = Actionability("civic", listOf(cancerType), drugs, level, significance,
                                                                "Predictive", HmfLevel(csvRecord["Level"]))

    companion object {
        private fun readLevel(levelField: String): String {
            return if (levelField.startsWith("R")) levelField.substring(1) else levelField
        }
    }
}
