package com.hartwig.hmftools.knowledgebaseimporter.oncoKb.readers

import com.hartwig.hmftools.knowledgebaseimporter.FusionReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoKbInput
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionEvent

object OncoFusionReader : SomaticEventReader<OncoKbInput, FusionEvent> {
    private val FUSION_SEPARATORS = setOf("-", " - ", "?")
    private val fusionReader = FusionReader(separators = FUSION_SEPARATORS)

    private fun matches(event: OncoKbInput) = event.variant.contains("Fusion".toRegex())

    override fun read(event: OncoKbInput): List<FusionEvent> {
        if (matches(event)) return listOfNotNull(fusionReader.read(event.gene, event.variant))
        return emptyList()
    }
}
