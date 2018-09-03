package com.hartwig.hmftools.knowledgebaseimporter.civic.readers

import com.hartwig.hmftools.knowledgebaseimporter.FusionReader
import com.hartwig.hmftools.knowledgebaseimporter.civic.input.CivicVariantInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionEvent

object CivicFusionReader : SomaticEventReader<CivicVariantInput, FusionEvent> {
    private val fusionReader = FusionReader(separators = setOf("-"))

    private fun matches(event: CivicVariantInput) = event.variantTypes.isNotEmpty() && event.variantTypes.all { it.contains("fusion") }

    override fun read(event: CivicVariantInput): List<FusionEvent> {
        if (matches(event)) return listOfNotNull(fusionReader.read(event.gene, event.variant))
        return emptyList()
    }
}
