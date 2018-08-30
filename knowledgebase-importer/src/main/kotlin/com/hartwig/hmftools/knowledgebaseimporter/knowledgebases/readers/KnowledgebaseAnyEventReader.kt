package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.KnowledgebaseEvent

// MIVO: reader interface for knowledgebase events that contain multiple variants, any of which would lead to actionability
// e.g:  V600E/V600K; V600E or V600K, etc
interface KnowledgebaseAnyEventReader<R, out T> : SomaticEventReader<R, T> where R : KnowledgebaseEvent, T : SomaticEvent {
    // MIVO: function that should split the single event into multiple events
    // e.g:  V600E/V600K -> [V600E, V600K]
    fun mapper(input: R): List<R>

    val nestedEventReader: KnowledgebaseEventReader<R, T>

    override fun read(event: R): List<T> {
        val nestedEvents = mapper(event)
        val somaticEvents = nestedEvents.map { nestedEvent -> nestedEventReader.read(nestedEvent) }
        return somaticEvents.flatten()
    }
}
