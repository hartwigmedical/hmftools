package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnowledgebaseEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import org.apache.logging.log4j.LogManager

data class KnowledgebaseEventReader<R : KnowledgebaseEvent>(val readers: List<SomaticEventReader<R, *>>) {
    private val logger = LogManager.getLogger("KnowledgebaseEventReader")

    fun read(event: R): List<SomaticEvent> {
        val events = readers.map { Pair(it.javaClass.name, it.read(event)) }.filterNot { it.second.isEmpty() }
        if (events.size > 1) {
            logger.warn("More than 1 reader (${events.joinToString(", ")}) returned events for record: $event")
        }
        if (events.isEmpty()) {
            logger.warn("Could not extract somatic event from:\t${event.source}\t${event.gene}\t${event.variant}\t${event.types}")
        }
        return events.map { it.second }.flatten()
    }
}
