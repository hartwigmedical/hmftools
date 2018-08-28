package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.GeneMutations
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader

object IclusionGeneMutationReader : SomaticEventReader<IclusionEvent, GeneMutations> {
    private fun match(event: IclusionEvent) = event.variant == "MUTATION" || event.variant == "ALTERATION"

    override fun read(event: IclusionEvent): List<GeneMutations> {
        if (match(event)) return listOf(GeneMutations(event.gene, event.transcript))
        return emptyList()
    }
}
