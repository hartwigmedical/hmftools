package com.hartwig.hmftools.knowledgebaseimporter.civic.readers

import com.hartwig.hmftools.knowledgebaseimporter.civic.input.CivicVariantInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.GenericRangeMutations
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers.CodonMatcher

object CivicRangeMutationReader : SomaticEventReader<CivicVariantInput, GenericRangeMutations> {
    private val MUTATION = "Mutation".toRegex(RegexOption.IGNORE_CASE)
    private val EXON_MUTATION = "Exon\\s[\\d]+(?:-[\\d]+)?\\sMutation".toRegex(RegexOption.IGNORE_CASE)
    private val DOMAIN_MUTATION = ".+\\sDomain\\sMutation".toRegex(RegexOption.IGNORE_CASE)
    private val G12G13 = "G12/G13".toRegex(RegexOption.IGNORE_CASE)
    private val variantPatterns = listOf(MUTATION, EXON_MUTATION, DOMAIN_MUTATION, G12G13)

    private fun match(event: CivicVariantInput) = event.hasPosition && !event.hasRefOrAlt &&
            (CodonMatcher.matches(event.variant) || variantPatterns.any { event.variant.matches(it) })

    override fun read(event: CivicVariantInput): List<GenericRangeMutations> {
        if (match(event)) {
            return listOf(GenericRangeMutations(event.gene, event.transcript.substringBefore("."), event.start.toInt(),
                                                event.stop.toInt()))
        }
        return emptyList()
    }
}
