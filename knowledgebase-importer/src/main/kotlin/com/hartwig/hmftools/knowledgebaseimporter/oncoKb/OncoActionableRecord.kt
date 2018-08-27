package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordMetadata
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfDrug
import org.apache.logging.log4j.LogManager

data class OncoActionableRecord(private val metadata: RecordMetadata, override val events: List<SomaticEvent>,
                                override val actionability: List<Actionability>) : RecordMetadata by metadata, ActionableRecord {
    companion object {
        private val logger = LogManager.getLogger("OncoActionableRecord")
        private val somaticEventReader = OncoSomaticEventReader()

        operator fun invoke(input: OncoActionableInput, treatmentTypeMap: Map<String, String>): OncoActionableRecord {
            val drugs = readDrugEntries(input, treatmentTypeMap)
            val actionability = Actionability("oncoKb", input.reference, listOf(input.`Cancer Type`), drugs, input.level,
                                              input.significance.name, "Predictive", input.hmfLevel, input.significance)
            val metadata = OncoMetadata(input.Gene, input.transcript)
            val events = somaticEventReader.read(input.Gene, input.transcript, input.Alteration)
            if (events.isEmpty()) {
                val aOrBLevelCount = actionability.filter { it.hmfLevel == "A" || it.hmfLevel == "B" }.size
                logger.warn("Could not extract somatic event from:\toncoKb\t${input.Gene}\t${input.Alteration}\t\t$aOrBLevelCount")
            }
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
