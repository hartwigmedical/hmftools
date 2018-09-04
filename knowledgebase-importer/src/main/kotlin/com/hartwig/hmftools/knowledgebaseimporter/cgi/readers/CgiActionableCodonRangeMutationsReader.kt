package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CodonRangeMutations
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader

object CgiActionableCodonRangeMutationsReader : SomaticEventReader<CgiActionableInput, CodonRangeMutations> {
    private val RANGE_PATTERN = "[0-9]+-[0-9]+".toRegex()

    private fun matches(event: CgiActionableInput) = event.`Alteration type` == "MUT" && RANGE_PATTERN.matches(event.variant)

    override fun read(event: CgiActionableInput): List<CodonRangeMutations> {
        if (matches(event)) {
            val range = codonRange(event)
            return listOf(CodonRangeMutations(event.gene, event.transcript, range.start, range.endInclusive))
        }
        return emptyList()
    }

    private fun codonRange(event: CgiActionableInput) =
            IntRange(event.variant.substringBefore("-").toInt(), event.variant.substringAfter("-").toInt())
}
