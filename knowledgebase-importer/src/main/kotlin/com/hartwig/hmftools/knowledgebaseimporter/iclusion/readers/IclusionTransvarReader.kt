package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.HgvsAnnotation
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader
import com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers.Matcher
import com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers.TransvarMatcher

object IclusionTransvarReader : SomaticEventReader<IclusionEvent, HgvsAnnotation> {
    private val nestedReaders = listOf(IclusionCDnaAnnotationReader, IclusionProteinAnnotationReader)

    private fun match(event: IclusionEvent) = TransvarMatcher.matches(event.variant)

    private fun contains(event: IclusionEvent) = TransvarMatcher.contains(event.variant)

    override fun read(event: IclusionEvent): List<HgvsAnnotation> {
        if (match(event)) return nestedReaders.flatMap { it.read(event) }
        if (contains(event)) {
            val events = event.variant.split(Matcher.SPLIT_PATTERN.toRegex()).map { event.copy(variant = it) }
            return events.flatMap { subEvent -> nestedReaders.flatMap { it.read(subEvent) } }
        }
        return emptyList()
    }
}
