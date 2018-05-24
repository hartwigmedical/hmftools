package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ActionableRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordMetadata
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.Actionability
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfLevel
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfResponse
import org.apache.commons.csv.CSVRecord

data class OncoActionableRecord(private val metadata: RecordMetadata, override val events: List<SomaticEvent>,
                                override val actionability: List<Actionability>) : RecordMetadata by metadata, ActionableRecord {
    companion object {
        private val somaticEventReader = OncoSomaticEventReader()

        operator fun invoke(record: CSVRecord): OncoActionableRecord {
            val gene = record["Gene"]
            val transcript = record["Isoform"]
            val level: String = readLevel(record["Level"])
            val significance = if (record["Level"].startsWith("R")) HmfResponse.Resistant else HmfResponse.Responsive
            val drugs = record["Drugs(s)"].split(",").map { it.trim() }
            val cancerType: String = record["Cancer Type"]
            val actionability = Actionability("oncoKb", listOf(cancerType), drugs, level, significance.name,
                                              "Predictive", HmfLevel(record["Level"]), significance)
            val alteration = record["Alteration"]
            val metadata = OncoMetadata(gene, transcript)
            return OncoActionableRecord(metadata, somaticEventReader.read(gene, transcript, alteration), actionability)
        }

        private fun readLevel(levelField: String): String = if (levelField.startsWith("R")) levelField.drop(1) else levelField
    }
}
