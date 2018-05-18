package com.hartwig.hmftools.knowledgebaseimporter.civic

import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfLevel
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfResponse
import org.apache.commons.csv.CSVRecord

data class CivicEvidence(private val csvRecord: CSVRecord, private val drugInteractionMap: Map<Int, String>) {
    val cancerType: String = csvRecord["disease"].orEmpty()
    val drugInteractionType = getDrugInteraction(csvRecord["evidence_id"].toIntOrNull(), drugInteractionMap)
    val doid: String = csvRecord["doid"].orEmpty()
    val drugs: List<String> = getDrugs(csvRecord["drugs"].orEmpty(), drugInteractionType)
    val type: String = csvRecord["evidence_type"].orEmpty()
    val direction: String = csvRecord["evidence_direction"].orEmpty()
    val level: String = csvRecord["evidence_level"].orEmpty()
    val significance: String = csvRecord["clinical_significance"].orEmpty()
    val actionabilityItems: List<Actionability> = Actionability("civic", listOf(cancerType), drugs, level, significance, type,
                                                                HmfLevel(level), HmfResponse(significance))

    companion object {
        private fun getDrugInteraction(evidenceId: Int?, drugInteractionMap: Map<Int, String>): String {
            evidenceId ?: return ""
            return drugInteractionMap[evidenceId].orEmpty()
        }

        private fun getDrugs(drugField: String, drugInteraction: String): List<String> {
            val drugs = drugField.split(",").map { it.trim() }.filterNot { it.isEmpty() }
            return when {
                drugInteraction.isBlank()                      -> drugs
                drugInteraction.toLowerCase() == "substitutes" -> drugs
                drugInteraction.toLowerCase() == "combination" -> listOf(drugs.joinToString(" + "))
                else                                           -> listOf(drugs.joinToString(", ") + " ($drugInteraction)")
            }
        }
    }
}
