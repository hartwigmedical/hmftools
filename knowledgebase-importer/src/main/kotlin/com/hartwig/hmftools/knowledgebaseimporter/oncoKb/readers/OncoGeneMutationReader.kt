package com.hartwig.hmftools.knowledgebaseimporter.oncoKb.readers

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.GeneMutations
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoKbInput

object OncoGeneMutationReader : SomaticEventReader<OncoKbInput, GeneMutations> {
    override fun read(event: OncoKbInput): List<GeneMutations> {
        return when (event.variant) {
            "Oncogenic Mutations" -> listOf(GeneMutations(event.gene, event.transcript))
            else                  -> emptyList()
        }
    }
}
