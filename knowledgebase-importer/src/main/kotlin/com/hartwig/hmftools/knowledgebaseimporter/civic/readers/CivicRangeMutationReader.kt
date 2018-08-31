package com.hartwig.hmftools.knowledgebaseimporter.civic.readers

import com.hartwig.hmftools.knowledgebaseimporter.civic.input.CivicVariantInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.GenericRangeMutations
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader

object CivicRangeMutationReader : SomaticEventReader<CivicVariantInput, GenericRangeMutations> {
    private const val RANGE_VARIANTS = "gene_variant|transcript_variant|exon_variant|coding_sequence_variant|protein_altering_variant"
    private val AND_PATTERN = "(\\+|and)".toRegex(RegexOption.IGNORE_CASE)

    private fun match(event: CivicVariantInput) = event.hasPosition && !event.hasRefOrAlt &&
            (isGenericMutation(event) || isGenericMissense(event))

    private fun isGenericMutation(event: CivicVariantInput) = event.variant.toLowerCase() == "mutation" ||
            event.variantTypes.any { it.contains(RANGE_VARIANTS.toRegex()) }

    private fun isGenericMissense(event: CivicVariantInput) = !event.variant.contains(AND_PATTERN) && event.variantTypes.size == 1
            && event.variantTypes.first() == "missense_variant"

    override fun read(event: CivicVariantInput): List<GenericRangeMutations> {
        if (match(event)) return listOf(GenericRangeMutations(event.gene, event.transcript.substringBefore("."), event.start.toInt(),
                                                              event.stop.toInt()))
        return emptyList()
    }
}
