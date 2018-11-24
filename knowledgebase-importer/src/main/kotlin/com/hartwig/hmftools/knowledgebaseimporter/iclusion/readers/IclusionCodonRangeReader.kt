package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CodonRangeMutations
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader

object IclusionCodonRangeReader : SomaticEventReader<IclusionEvent, CodonRangeMutations> {
    private val EXON_RANGE_MUTATION_PATTERN = "G12/G13".toRegex()

    private fun matches(event: IclusionEvent) = IclusionCodonRangeReader.EXON_RANGE_MUTATION_PATTERN.matches(event.variant)

    override fun read(event: IclusionEvent): List<CodonRangeMutations> {
        if (IclusionCodonRangeReader.matches(event)) {
            val range = codonRange(event)
            return listOf(CodonRangeMutations(event.gene, event.transcript, range.start, range.endInclusive))
        }
        return emptyList()
    }

    private fun codonRange(event: IclusionEvent) =
            // TODO (LISC): Can probably go, pending confirmation from Paul.
            IntRange(12, 13)

}