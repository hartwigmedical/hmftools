package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnownRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordMetadata
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import org.apache.commons.csv.CSVRecord

data class OncoKnownRecord(private val metadata: RecordMetadata, override val additionalInfo: String,
                           override val events: List<SomaticEvent>) : KnownRecord, RecordMetadata by metadata {
    companion object {
        private val somaticEventReader = OncoSomaticEventReader()

        operator fun invoke(record: CSVRecord): OncoKnownRecord {
            val (gene, alteration) = correctRecord(record["Gene"], record["Alteration"])
            val transcript = record["Isoform"]
            val metadata = OncoMetadata(gene, transcript)
            return OncoKnownRecord(metadata, record["Oncogenicity"], somaticEventReader.read(gene, transcript, alteration))
        }

        private fun correctRecord(gene: String, alteration: String): Pair<String, String> = when {
            alteration.contains(Regex("IGH-NKX2")) && gene == "NKX2-1" ->
                Pair(gene, alteration.replace("IGH-NKX2", "IGH-NKX2-1"))
            else                                                       -> Pair(gene, alteration)
        }
    }
}
