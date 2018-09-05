package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader

object CgiActionableVariantReader : SomaticEventReader<CgiActionableInput, SomaticEvent> {
    private val readers = listOf(CgiGDnaReader, CgiActionableCDnaAnnotationReader, CgiActionableProteinAnnotationReader,
                                 CgiActionableCodonMutationsReader, CgiActionableCodonRangeMutationsReader,
                                 CgiActionableGeneMutationsReader, CgiActionableMutationListReader)

    private fun match(event: CgiActionableInput) = event.`Alteration type` == "MUT"

    override fun read(event: CgiActionableInput): List<SomaticEvent> {
        if (match(event)) return readers.flatMap { it.read(event) }
        return emptyList()
    }
}
