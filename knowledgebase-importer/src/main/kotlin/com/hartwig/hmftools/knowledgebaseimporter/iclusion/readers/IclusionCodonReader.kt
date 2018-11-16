package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CodonMutations
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers.CodonMatcher

object IclusionCodonReader : SomaticEventReader<IclusionEvent, CodonMutations> {
    private const val EXON_RANGE_MUTATION_PATTERN = "G12+\\s*-\\s*G13"

    override fun read(event: IclusionEvent): List<CodonMutations> {
        return when {
            CodonMatcher.matches(event.variant) -> { listOf(CodonMutations(event.gene, event.transcript, codonNumber(event)))}
            CodonMatcher.matches(EXON_RANGE_MUTATION_PATTERN) -> { listOf(CodonMutations(event.gene, event.transcript, codonNumber(event)))}
            else -> emptyList()
        }
    }

    private fun codonNumber(event: IclusionEvent) = "([\\d]+)".toRegex().find(event.variant)!!.groupValues[1].toInt()
}
