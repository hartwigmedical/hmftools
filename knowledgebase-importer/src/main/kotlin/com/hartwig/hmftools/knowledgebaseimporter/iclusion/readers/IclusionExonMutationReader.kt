package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.EventType
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ExonMutations
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader

object IclusionExonMutationReader : SomaticEventReader<IclusionEvent, ExonMutations> {
    private const val SINGLE_EXON_MUTATION_PATTERN = "exon\\s+[0-9]+\\s+mutation"
    private const val EXON_RANGE_MUTATION_PATTERN = "exon\\s+[0-9]+\\s*-\\s*[0-9]+\\s+mutation"
    private val SINGLE_EXON_REGEX = SINGLE_EXON_MUTATION_PATTERN.toRegex(RegexOption.IGNORE_CASE)
    private val EXON_RANGE_REGEX = EXON_RANGE_MUTATION_PATTERN.toRegex(RegexOption.IGNORE_CASE)

    val EXON_MUTATION_REGEX = "(?:$SINGLE_EXON_MUTATION_PATTERN|$EXON_RANGE_MUTATION_PATTERN)".toRegex(RegexOption.IGNORE_CASE)

    private fun match(event: IclusionEvent): Boolean {
        return event.types.size == 1 && event.types.contains(EventType.MUT) && event.variant.matches(EXON_MUTATION_REGEX)
    }

    override fun read(event: IclusionEvent): List<ExonMutations> {
        if (match(event)) return readMutations(event)
        return emptyList()
    }

    private fun readMutations(event: IclusionEvent): List<ExonMutations> {
        return when {
            event.variant.matches(SINGLE_EXON_REGEX) -> {
                val exonNumber = "[0-9]+".toRegex().find(event.variant)!!.value
                listOf(ExonMutations(event.gene, event.transcript, exonNumber.toInt()))
            }
            event.variant.matches(EXON_RANGE_REGEX)  -> {
                val exonNumbers = "([0-9]+)".toRegex().findAll(event.variant).toList().map { it.value }
                val exonRange = IntRange(exonNumbers[0].toInt(), exonNumbers[1].toInt())
                exonRange.map { exonNumber -> ExonMutations(event.gene, event.transcript, exonNumber) }
            }
            else                                     -> emptyList()
        }
    }
}
