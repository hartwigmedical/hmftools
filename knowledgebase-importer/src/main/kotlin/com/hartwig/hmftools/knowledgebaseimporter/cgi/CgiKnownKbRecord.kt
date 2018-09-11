package com.hartwig.hmftools.knowledgebaseimporter.cgi

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiKnownInput
import com.hartwig.hmftools.knowledgebaseimporter.cgi.readers.CgiKnownInputReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnownRecord
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.RecordMetadata
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.KnowledgebaseEventReader

data class CgiKnownKbRecord(private val metadata: RecordMetadata, override val additionalInfo: String,
                            override val events: List<SomaticEvent>) : RecordMetadata by metadata, KnownRecord {
    companion object {
        private val reader = KnowledgebaseEventReader("cgi", CgiKnownInputReader)

        operator fun invoke(input: CgiKnownInput): CgiKnownKbRecord {
            val metadata = CgiMetadata(input.gene, input.transcript)
            val events = reader.read(input)
            return CgiKnownKbRecord(metadata, "", events)
        }
    }
}
