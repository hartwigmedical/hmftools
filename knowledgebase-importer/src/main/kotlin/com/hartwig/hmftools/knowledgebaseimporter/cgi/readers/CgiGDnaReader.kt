package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.GDnaVariant
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader

object CgiGDnaReader : SomaticEventReader<CgiInput, GDnaVariant> {
    override fun read(event: CgiInput): List<GDnaVariant> =
            event.gdna.split("__").map { it.trim() }.filterNot { it.isBlank() }.map { GDnaVariant(it) }
}
