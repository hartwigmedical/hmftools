package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CDnaAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.CDnaAnnotationReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader

object CgiActionableCDnaAnnotationReader : SomaticEventReader<CgiActionableInput, CDnaAnnotation> {
    override fun read(event: CgiActionableInput): List<CDnaAnnotation> {
        if (event.`Alteration type` == "MUT")
            return CDnaAnnotationReader.read(event.copy(Alteration = event.cDNA.orEmpty()))
        return emptyList()
    }
}
