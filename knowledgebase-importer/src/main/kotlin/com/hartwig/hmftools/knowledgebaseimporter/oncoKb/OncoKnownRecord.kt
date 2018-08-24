package com.hartwig.hmftools.knowledgebaseimporter.oncoKb

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnownRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordMetadata
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoKnownInput

data class OncoKnownRecord(private val metadata: RecordMetadata, override val additionalInfo: String,
                           override val events: List<SomaticEvent>) : KnownRecord, RecordMetadata by metadata {
    companion object {
        private val somaticEventReader = OncoSomaticEventReader()

        operator fun invoke(input: OncoKnownInput): OncoKnownRecord {
            val metadata = OncoMetadata(input.Gene, input.transcript)
            return OncoKnownRecord(metadata, input.Oncogenicity, somaticEventReader.read(input.Gene, input.transcript, input.Alteration,
                                                                                         input.`Mutation Effect`))
        }
    }
}
