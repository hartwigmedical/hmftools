package com.hartwig.hmftools.knowledgebaseimporter.civic.readers

import com.hartwig.hmftools.knowledgebaseimporter.FusionReader
import com.hartwig.hmftools.knowledgebaseimporter.civic.input.CivicVariantInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.OtherEvents
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader

object CivicMultipleEventsReader : SomaticEventReader<CivicVariantInput, SomaticEvent> {
    private val fusionReader = FusionReader(separators = setOf("-"))

    override fun read(event: CivicVariantInput): List<OtherEvents> {
        val events = mutableMapOf<String, List<SomaticEvent>>()
        if (containsFusion(event)) events["fusion"] = listOfNotNull(fusionReader.read(event.gene, event.variant))
        val variants = CivicVariantReader.read(event.copy(variant_types = ""))
        if (variants.isNotEmpty()) events["variants"] = variants
        return buildEvents(events)
    }

    private fun containsFusion(event: CivicVariantInput) =
            event.variantTypes.any { it.contains("fusion") } || event.variant.contains("fusion".toRegex(RegexOption.IGNORE_CASE))

    private fun buildEvents(eventMap: Map<String, List<SomaticEvent>>): List<OtherEvents> {
        return if (eventMap.filterValues { it.isNotEmpty() }.size < 2) emptyList()
        else emptyList<List<SomaticEvent>>().combine(eventMap).map { OtherEvents(it) }
    }

    private fun <T> List<List<T>>.combine(map: Map<String, List<T>>): List<List<T>> {
        if (map.isEmpty()) return this
        val key = map.keys.first()
        return map[key]!!.map { value ->
            if (this.isEmpty()) listOf(value)
            else this.flatMap { it + value }
        }.combine(map.filterKeys { it != key })
    }
}
