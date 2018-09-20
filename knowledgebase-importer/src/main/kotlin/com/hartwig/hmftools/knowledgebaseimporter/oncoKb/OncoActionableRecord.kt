package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordMetadata
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.KnowledgebaseEventReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.ProteinAnnotationReader
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.readers.*
import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfDrug

data class OncoActionableRecord(private val metadata: RecordMetadata, override val events: List<SomaticEvent>,
                                override val actionability: List<Actionability>) :
        RecordMetadata by metadata, ActionableRecord {
    companion object {
        private val multiSnvReader = OncoAnyActionableEventReader(KnowledgebaseEventReader("oncoKb", ProteinAnnotationReader))
        private val reader = KnowledgebaseEventReader("oncoKb", OncoCnvReader, OncoExonMutationReader, OncoFusionReader,
                                                      OncoGeneMutationReader, ProteinAnnotationReader, multiSnvReader)

        operator fun invoke(input: OncoActionableInput, treatmentTypeMap: Map<String, String>): OncoActionableRecord {
            val drugs = readDrugEntries(input, treatmentTypeMap)
            val actionability = Actionability("oncoKb", input.reference, listOf(input.`Cancer Type`), drugs, input.level,
                                              input.hmfResponse.name, "Predictive", input.hmfLevel, input.hmfResponse)
            val metadata = OncoMetadata(input.gene, input.transcript)
            val events = reader.read(input)
            return OncoActionableRecord(metadata, events, actionability)
        }

        private fun readDrugEntries(input: OncoActionableInput, treatmentTypeMap: Map<String, String>): List<HmfDrug> {
            val drugNames = input.`Drugs(s)`.split(",").map { it.trim() }
            return drugNames.map { annotateDrugEntry(it, treatmentTypeMap) }
        }

        private fun annotateDrugEntry(entry: String, treatmentTypeMap: Map<String, String>): HmfDrug {
            val entryType = entry.split("+")
                    .joinToString(" + ") { treatmentTypeMap[it.trim().toLowerCase()] ?: "Unknown" }
            return HmfDrug(entry, entryType)
        }
    }
}
