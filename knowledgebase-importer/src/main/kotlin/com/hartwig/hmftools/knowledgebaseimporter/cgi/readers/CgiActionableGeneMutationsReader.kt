package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.GeneMutations
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader

object CgiActionableGeneMutationsReader : SomaticEventReader<CgiActionableInput, GeneMutations> {
    override fun read(event: CgiActionableInput): List<GeneMutations> {
        if (event.`Alteration type` == "MUT" && event.variant == ".") return listOf(GeneMutations(event.gene, event.transcript))
        return emptyList()
    }
}
