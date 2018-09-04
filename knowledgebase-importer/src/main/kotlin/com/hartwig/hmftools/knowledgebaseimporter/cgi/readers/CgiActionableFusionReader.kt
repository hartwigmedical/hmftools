package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.FusionReader
import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionEvent

object CgiActionableFusionReader : SomaticEventReader<CgiActionableInput, FusionEvent> {
    private val FUSION_SEPARATORS = setOf("__")
    private val fusionReader = FusionReader(separators = FUSION_SEPARATORS)

    override fun read(event: CgiActionableInput): List<FusionEvent> {
        if (event.`Alteration type` == "FUS") return listOfNotNull(fusionReader.read(event.gene, event.Alteration))
        return emptyList()
    }
}
