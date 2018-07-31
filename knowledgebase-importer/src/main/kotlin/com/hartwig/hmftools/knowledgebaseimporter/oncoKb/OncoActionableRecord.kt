package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordMetadata
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfDrug
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfLevel
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfResponse
import org.apache.commons.csv.CSVRecord
import org.apache.logging.log4j.LogManager

data class OncoActionableRecord(private val metadata: RecordMetadata, override val events: List<SomaticEvent>,
                                override val actionability: List<Actionability>) : RecordMetadata by metadata, ActionableRecord {
    companion object {
        private val logger = LogManager.getLogger("OncoActionableRecord")
        private val somaticEventReader = OncoSomaticEventReader()

        operator fun invoke(record: CSVRecord, treatmentTypeMap: Map<String, String>): OncoActionableRecord {
            val gene = record["Gene"]
            val transcript = record["Isoform"]
            val level: String = readLevel(record["Level"])
            val significance = if (record["Level"].startsWith("R")) HmfResponse.Resistant else HmfResponse.Responsive
            val drugs = readDrugEntries(record, treatmentTypeMap)
            val cancerType: String = readCancerType(record)
            val alteration = record["Alteration"]
            val actionability = Actionability("oncoKb", "$gene $alteration", listOf(cancerType), drugs, level,
                                              significance.name, "Predictive", HmfLevel(record["Level"]), significance)
            val metadata = OncoMetadata(gene, transcript)
            val events = somaticEventReader.read(gene, transcript, alteration)
            if (events.isEmpty()) {
                val aOrBLevelCount = actionability.filter { it.hmfLevel == "A" || it.hmfLevel == "B" }.size
                logger.warn("Could not extract somatic event from:\toncoKb\t$gene\t$alteration\t\t$aOrBLevelCount")
            }
            return OncoActionableRecord(metadata, events, actionability)
        }

        private fun readLevel(levelField: String): String = if (levelField.startsWith("R")) levelField.drop(1) else levelField

        private fun readDrugEntries(record: CSVRecord, treatmentTypeMap: Map<String, String>): List<HmfDrug> {
            val drugNames = record["Drugs(s)"].split(",").map { it.trim() }
            return drugNames.map { annotateDrugEntry(it, treatmentTypeMap) }
        }

        private fun annotateDrugEntry(entry: String, treatmentTypeMap: Map<String, String>): HmfDrug {
            val entryType = entry.split("+")
                    .joinToString(" + ") { treatmentTypeMap[it.trim().toLowerCase()] ?: "Unknown" }
            return HmfDrug(entry, entryType)
        }

        private fun readCancerType(record: CSVRecord): String {
            val cancerType = record["Cancer Type"]
            return if (cancerType == "Melanoma") "Skin Melanoma" else cancerType
        }
    }
}
