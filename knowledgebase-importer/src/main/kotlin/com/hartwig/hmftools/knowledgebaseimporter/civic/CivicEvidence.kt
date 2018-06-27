package com.hartwig.hmftools.knowledgebaseimporter.civic

import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfDrug
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfLevel
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfResponse
import org.apache.commons.csv.CSVRecord

data class CivicEvidence(private val csvRecord: CSVRecord, private val drugInteractionMap: Map<Int, String>,
                         private val treatmentTypeMap: Map<String, String>) {
    val cancerType: String = csvRecord["disease"].orEmpty()
    val drugInteractionType = getDrugInteraction(csvRecord["evidence_id"].toIntOrNull(), drugInteractionMap)
    val doid: String = csvRecord["doid"].orEmpty()
    val drugs: List<HmfDrug> = getDrugs(csvRecord["drugs"].orEmpty(), drugInteractionType, treatmentTypeMap)
    val type: String = csvRecord["evidence_type"].orEmpty()
    val direction: String = csvRecord["evidence_direction"].orEmpty()
    val level: String = csvRecord["evidence_level"].orEmpty()
    val significance: String = csvRecord["clinical_significance"].orEmpty()
    val actionabilityItems: List<Actionability> = Actionability("civic", csvRecord["evidence_id"], listOf(cancerType), drugs, level,
                                                                significance, type, HmfLevel(level), HmfResponse(significance))

    companion object {
        private const val drugCharacters = "A-Z0-9a-z-\\s"
        private const val squareBracketPart = "\\[.*?\\]"
        private const val roundBracketPart = "\\(.*?\\)"

        private fun getDrugInteraction(evidenceId: Int?, drugInteractionMap: Map<Int, String>): String {
            evidenceId ?: return ""
            return drugInteractionMap[evidenceId].orEmpty()
        }

        private fun getDrugs(drugField: String, drugInteraction: String, treatmentTypeMap: Map<String, String>): List<HmfDrug> {
            val drugs = splitDrugs(drugField).filterNot { it.isEmpty() }
            val treatmentTypes = drugs.map { treatmentTypeMap[it.toLowerCase()] ?: "Unknown" }
            return when {
                drugInteraction.isBlank()                      -> drugs.zip(treatmentTypes).map { HmfDrug(it.first, it.second) }
                drugInteraction.toLowerCase() == "substitutes" -> drugs.zip(treatmentTypes).map { HmfDrug(it.first, it.second) }
                drugInteraction.toLowerCase() == "combination" -> {
                    val name = drugs.joinToString(" + ")
                    val type = treatmentTypes.joinToString(" + ")
                    listOf(HmfDrug(name, type))
                }
                else                                           -> {
                    val name = drugs.joinToString(", ") + " ($drugInteraction)"
                    val type = treatmentTypes.joinToString(", ") + " ($drugInteraction)"
                    listOf(HmfDrug(name, type))
                }
            }
        }

        fun splitDrugs(drugField: String): List<String> {
            val drugNamePattern = "([$drugCharacters]+(?:\\s*(?:$roundBracketPart|$squareBracketPart)[$drugCharacters]*)*)"
            val matchResult = drugNamePattern.toRegex().findAll(drugField)
            return matchResult.map { it.groupValues[1].trim() }.toList()
        }
    }
}
