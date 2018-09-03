package com.hartwig.hmftools.knowledgebaseimporter.civic.readers

import com.hartwig.hmftools.knowledgebaseimporter.civic.input.CivicVariantInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader

object CivicVariantReader : SomaticEventReader<CivicVariantInput, SomaticEvent> {
    private val nestedReaders = listOf(CivicKnowledgebaseVariantReader, CivicCDnaAnnotationReader)

    fun match(event: CivicVariantInput) = event.variantTypes.none { it.contains("fusion") }

    override fun read(event: CivicVariantInput): List<SomaticEvent> {
        return if (match(event)) nestedReaders.flatMap { it.read(event) }
        else emptyList()
    }
}
