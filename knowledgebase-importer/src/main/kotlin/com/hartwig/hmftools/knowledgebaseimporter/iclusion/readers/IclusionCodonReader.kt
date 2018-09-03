package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CodonMutations
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers.CodonMatcher

object IclusionCodonReader : SomaticEventReader<IclusionEvent, CodonMutations> {
    override fun read(event: IclusionEvent): List<CodonMutations> {
        return if (CodonMatcher.matches(event.variant)) listOf(CodonMutations(event.gene, event.transcript, codonNumber(event)))
        else emptyList()
    }

    private fun codonNumber(event: IclusionEvent): Int {
        return "([\\d]+)".toRegex().find(event.variant)!!.groupValues[1].toInt()
    }
}
