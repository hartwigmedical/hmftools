package com.hartwig.hmftools.knowledgebaseimporter.civic.readers

import com.hartwig.hmftools.knowledgebaseimporter.civic.input.CivicVariantInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnowledgebaseVariant
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader

object CivicKnowledgebaseVariantReader : SomaticEventReader<CivicVariantInput, KnowledgebaseVariant> {
    // MIVO: allow ref or alt to be empty since civic annotates insertions as ['' -> 'GATC'] and deletions as ['GATC' -> '']
    private fun match(event: CivicVariantInput) = event.hasPosition && event.hasRefOrAlt

    override fun read(event: CivicVariantInput): List<KnowledgebaseVariant> {
        return if (match(event)) {
            val ref = if (event.reference_bases.isBlank()) null else event.reference_bases
            val alt = if (event.variant_bases.isBlank()) null else event.variant_bases
            listOf(KnowledgebaseVariant(event.gene, event.chromosome, event.start.toLong(), ref, alt))
        } else {
            emptyList()
        }
    }
}
