package com.hartwig.hmftools.knowledgebaseimporter.cgi

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiKnownInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.*
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.SequenceVariantType

data class CgiKnownKbRecord(private val metadata: RecordMetadata, override val additionalInfo: String,
                            override val events: List<SomaticEvent>) : RecordMetadata by metadata, KnownRecord {
    companion object {
        operator fun invoke(input: CgiKnownInput): CgiKnownKbRecord {
            val metadata = CgiMetadata(input.gene, input.transcript)
            val events = readSomaticEvents(input)
            return CgiKnownKbRecord(metadata, "", events)
        }

        private fun readSomaticEvents(input: CgiKnownInput): List<SomaticEvent> {
            val gDnaVariants = input.gdna.split("__").map { it.trim() }.filterNot { it.isBlank() }.map { GDnaVariant(it) }
            val proteinAnnotations = listOfNotNull(input.protein).filterNot { it.isBlank() }.map {
                ProteinAnnotation(input.transcript, it, SequenceVariantType.OTHER)
            }
            return gDnaVariants + proteinAnnotations
        }
    }
}
