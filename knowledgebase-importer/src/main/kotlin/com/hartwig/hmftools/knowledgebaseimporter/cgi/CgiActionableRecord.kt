package com.hartwig.hmftools.knowledgebaseimporter.cgi

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.cgi.readers.CgiActionableCnvReader
import com.hartwig.hmftools.knowledgebaseimporter.cgi.readers.CgiActionableFusionReader
import com.hartwig.hmftools.knowledgebaseimporter.cgi.readers.CgiActionableVariantReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordMetadata
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.KnowledgebaseEventReader
import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfDrug
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfLevel
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfResponse

data class CgiActionableRecord(private val metadata: RecordMetadata, override val events: List<SomaticEvent>,
                               override val actionability: List<Actionability>, val cgiDrugs: List<CgiDrug>) : RecordMetadata by metadata,
        ActionableRecord {
    companion object {
        private val reader = KnowledgebaseEventReader("cgi", CgiActionableVariantReader, CgiActionableFusionReader, CgiActionableCnvReader)

        operator fun invoke(input: CgiActionableInput, treatmentTypeMap: Map<String, String>): CgiActionableRecord {
            val metadata = CgiMetadata(input.gene, input.transcript ?: "na")
            val events = reader.read(input)
            val actionability = readActionability(input, treatmentTypeMap)
            return CgiActionableRecord(metadata, events, actionability, readCgiDrugs(input))
        }

        private fun readActionability(input: CgiActionableInput, treatmentTypeMap: Map<String, String>): List<Actionability> {
            val cancerTypes = input.`Primary Tumor type`.split(";").map { it.trim() }
            val level = input.`Evidence level`
            val association = input.Association
            return Actionability("cgi", input.Alteration, cancerTypes, readDrugs(input, treatmentTypeMap), level, association,
                                 "Predictive", HmfLevel(level), HmfResponse(association))
        }

        private fun readDrugs(input: CgiActionableInput, treatmentTypeMap: Map<String, String>): List<HmfDrug> {
            val drugs = readDrugNames(input)
            return drugs.map { name ->
                if (name.contains("+")) {
                    val type = name.split("+").map { it.trim() }
                            .joinToString(" + ") { treatmentTypeMap[it.toLowerCase()] ?: "Unknown" }
                    HmfDrug(name, type)
                } else {
                    HmfDrug(name, treatmentTypeMap[name.toLowerCase()] ?: "Unknown")
                }
            }
        }

        private fun readDrugNames(input: CgiActionableInput): List<String> {
            val drugNames = readDrugsField(input.Drug.orEmpty())
            return if (drugNames.isEmpty()) {
                readDrugsField(input.`Drug family`.orEmpty())
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
    }
}
