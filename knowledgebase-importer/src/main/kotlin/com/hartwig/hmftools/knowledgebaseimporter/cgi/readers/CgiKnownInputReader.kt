package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiKnownInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.ProteinAnnotationReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader

object CgiKnownInputReader : SomaticEventReader<CgiKnownInput, SomaticEvent> {
    private val readers = listOf(CgiGDnaReader, ProteinAnnotationReader)
    override fun read(event: CgiKnownInput): List<SomaticEvent> = readers.flatMap { it.read(event) }
}
