package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.output.CnvEvent

object CgiActionableCnvReader : SomaticEventReader<CgiActionableInput, CnvEvent> {
    override fun read(event: CgiActionableInput): List<CnvEvent> {
        if (event.`Alteration type` == "CNA") {
            return if (event.Alteration.contains("amp")) listOf(CnvEvent.amplification(event.gene))
            else listOf(CnvEvent.deletion(event.gene))
        }
        return emptyList()
    }
}
