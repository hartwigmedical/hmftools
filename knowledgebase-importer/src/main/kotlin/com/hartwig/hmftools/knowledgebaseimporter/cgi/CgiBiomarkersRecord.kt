package com.hartwig.hmftools.knowledgebaseimporter.cgi

import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfLevel
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfResponse
import org.apache.commons.csv.CSVRecord

data class CgiBiomarkersRecord(private val csvRecord: CSVRecord) {
    val gene: String = csvRecord["Gene"]
    val gdna: String = csvRecord["gDNA"] ?: "na"
    val cdna: String = csvRecord["cDNA"] ?: "na"
    val alterationType: String = csvRecord["Alteration type"]
    val alteration: String = csvRecord["Alteration"]
    val transcript: String = csvRecord["transcript"] ?: "na"
    val protein: String = csvRecord["individual_mutation"]?.substringAfter(':', "") ?: "na"
    val drugs: List<String> = readDrugs(csvRecord)
    val level: String = csvRecord["Evidence level"]
    val association: String = csvRecord["Association"]
    val cancerTypes = csvRecord["Primary Tumor type"].split(";").map { it.trim() }
    val actionabilityItems: List<Actionability> = Actionability("cgi", cancerTypes, drugs, level, association, "Predictive",
                                                                highestLevel(level), HmfResponse(association))

    companion object {
        private fun readDrugs(csvRecord: CSVRecord): List<String> {
            val drugNames = readDrugsField(csvRecord["Drug"].orEmpty())
            return if (drugNames.isEmpty()) {
                readDrugsField(csvRecord["Drug family"].orEmpty())
            } else {
                drugNames
            }
        }

        private fun readDrugsField(drugField: String): List<String> {
            return drugField.replace(";", " + ")
                    .removeSurrounding("[", "]")
                    .split(",")
                    .filterNot { it.isBlank() }
        }

        private fun highestLevel(levelField: String): HmfLevel {
            return levelField.split(",").map { it.trim() }.map { HmfLevel(it) }.sorted().firstOrNull() ?: HmfLevel.UNKNOWN
        }
    }
}
