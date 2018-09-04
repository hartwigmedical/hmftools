package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ProteinAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.ProteinAnnotationReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader

object CgiActionableProteinAnnotationReader : SomaticEventReader<CgiActionableInput, ProteinAnnotation> {
    override fun read(event: CgiActionableInput): List<ProteinAnnotation> {
        if (event.`Alteration type` == "MUT")
            return ProteinAnnotationReader.read(event.copy(variant = event.individual_mutation.orEmpty().substringAfter(':', "")))
        return emptyList()
    }
}
