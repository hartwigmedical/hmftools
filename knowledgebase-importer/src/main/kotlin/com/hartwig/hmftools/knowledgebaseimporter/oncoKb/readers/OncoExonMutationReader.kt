package com.hartwig.hmftools.knowledgebaseimporter.oncoKb.readers

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ExonMutations
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoKbInput

object OncoExonMutationReader : SomaticEventReader<OncoKbInput, ExonMutations> {
    private val EXON_PATTERN = "Exon\\s?([0-9]+)\\s+Mutations".toRegex(RegexOption.IGNORE_CASE)

    override fun read(event: OncoKbInput): List<ExonMutations> {
        return when {
            event.variant.contains(EXON_PATTERN) -> listOf(ExonMutations(event.gene, event.transcript, extractExon(event.variant)))
            else                                 -> emptyList()
        }
    }

    private fun extractExon(alteration: String): Int {
        val matchResult = EXON_PATTERN.find(alteration)
        return matchResult!!.groupValues[1].toInt()
    }
}
