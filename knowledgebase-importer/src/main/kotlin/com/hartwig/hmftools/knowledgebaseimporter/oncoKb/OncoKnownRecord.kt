package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnownRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordMetadata
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.KnowledgebaseEventReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.ProteinAnnotationReader
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoKnownInput
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.readers.*
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionEvent

data class OncoKnownRecord(private val metadata: RecordMetadata, override val additionalInfo: String,
                           override val events: List<SomaticEvent>) : KnownRecord,
        RecordMetadata by metadata {
    companion object {
        private val multiSnvReader = OncoAnyKnownEventReader(KnowledgebaseEventReader("oncoKb", ProteinAnnotationReader))
        private val reader = KnowledgebaseEventReader("oncoKb", OncoCnvReader, OncoExonMutationReader, OncoFusionReader,
                                                      OncoGeneMutationReader, ProteinAnnotationReader, multiSnvReader)

        operator fun invoke(input: OncoKnownInput): OncoKnownRecord {
            val metadata = OncoMetadata(input.gene, input.transcript!!)
            val events = readEvents(input)
            return OncoKnownRecord(metadata, input.Oncogenicity, events)
        }

        private fun readEvents(input: OncoKnownInput): List<SomaticEvent> {
            val events = reader.read(input)
            return events.filterNot { it is FusionEvent && input.`Mutation Effect`.contains("Loss-of-function") }
        }
    }
}
