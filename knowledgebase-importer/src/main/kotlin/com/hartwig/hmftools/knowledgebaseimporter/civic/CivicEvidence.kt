package com.hartwig.hmftools.knowledgebaseimporter.civic

import com.hartwig.hmftools.knowledgebaseimporter.civic.input.CivicEvidenceInput
import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfDrug
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfLevel
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfResponse

data class CivicEvidence(private val input: CivicEvidenceInput, private val drugInteractionMap: Map<Int, String>,
                         private val treatmentTypeMap: Map<String, String>) {
    val disease: String = input.disease
    val doid: String = input.doid
    val type: String = input.evidence_type
    val direction: String = input.evidence_direction
    val level: String = input.evidence_level
    val actionabilityItems = Actionability("civic", input.evidence_id, listOf(disease),
                                           getDrugs(input, drugInteractionMap, treatmentTypeMap), level, input.clinical_significance, type,
                                           HmfLevel(level), HmfResponse(input.clinical_significance))

    companion object {
        private const val drugCharacters = "A-Z0-9a-z-\\s"
        private const val squareBracketPart = "\\[.*?\\]"
        private const val roundBracketPart = "\\(.*?\\)"

        private fun getDrugs(input: CivicEvidenceInput, drugInteractionMap: Map<Int, String>,
                             treatmentTypeMap: Map<String, String>): List<HmfDrug> {
            val drugInteraction = drugInteractionMap[input.evidence_id.toInt()].orEmpty()
            val drugNames = splitDrugs(input.drugs).filterNot { it.isEmpty() }
            return when {
                drugInteraction.isBlank()                      -> drugNames.map { toHmfDrug(it, treatmentTypeMap) }
                drugInteraction.toLowerCase() == "substitutes" -> drugNames.map { toHmfDrug(it, treatmentTypeMap) }
                drugInteraction.toLowerCase() == "combination" -> listOf(combiHmfDrug(drugNames, treatmentTypeMap))
                else                                           -> listOf(otherHmfDrug(drugNames, drugInteraction, treatmentTypeMap))
            }
        }

        fun splitDrugs(drugField: String): List<String> {
            val drugNamePattern = "([$drugCharacters]+(?:\\s*(?:$roundBracketPart|$squareBracketPart)[$drugCharacters]*)*)"
            val matchResult = drugNamePattern.toRegex().findAll(drugField)
            return matchResult.map { it.groupValues[1].trim() }.toList()
        }

        private fun drugType(name: String, typeMap: Map<String, String>): String = typeMap[name.toLowerCase()] ?: "Unknown"
        private fun toHmfDrug(name: String, typeMap: Map<String, String>): HmfDrug = HmfDrug(name, drugType(name, typeMap))

        private fun combiHmfDrug(names: List<String>, typeMap: Map<String, String>): HmfDrug {
            val types = names.map { drugType(it, typeMap) }
            return HmfDrug(names.joinToString(" + "), types.joinToString(" + "))
        }

        private fun otherHmfDrug(names: List<String>, interaction: String, typeMap: Map<String, String>): HmfDrug {
            val name = names.joinToString(", ") + " ($interaction)"
            val type = names.joinToString(", ") { drugType(it, typeMap) } + " ($interaction)"
            return HmfDrug(name, type)
        }
    }
}
