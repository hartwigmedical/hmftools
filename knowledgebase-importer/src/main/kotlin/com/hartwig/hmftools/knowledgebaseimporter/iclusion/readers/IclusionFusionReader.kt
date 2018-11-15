package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.FusionReader
import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.output.FusionEvent

object IclusionFusionReader : SomaticEventReader<IclusionEvent, FusionEvent> {
    private const val FUSION_SEPARATOR = "-"
    private val fusionReader = FusionReader(separators = setOf(FUSION_SEPARATOR))

    private fun match(event: IclusionEvent): Boolean {
        val variant = event.variant.toLowerCase()
        return variant == "fusions" || variant == "rearrangement" || variant.contains("fusion") ||
                (variant.contains(IclusionFusionReader.FUSION_SEPARATOR)
                && event.gene.length > 2 && variant.contains(event.gene.toLowerCase().substring(0, 3)))
    }

    override fun read(event: IclusionEvent): List<FusionEvent> {
        if (match(event)) return listOfNotNull(fusionReader.read(event.gene, event.variant))
        return emptyList()
    }
}
