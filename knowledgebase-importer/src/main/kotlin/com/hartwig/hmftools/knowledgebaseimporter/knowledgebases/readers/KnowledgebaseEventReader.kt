package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.KnowledgebaseEvent
import org.apache.logging.log4j.LogManager

data class KnowledgebaseEventReader<in R, out T>(val source: String, val readers: List<SomaticEventReader<R, T>>) : SomaticEventReader<R, T>
        where R : KnowledgebaseEvent, T : SomaticEvent {

    companion object {
        private val logger = LogManager.getLogger("KnowledgebaseEventReader")
        operator fun <R, T> invoke(source: String, vararg readers: SomaticEventReader<R, T>): KnowledgebaseEventReader<R, T>
                where R : KnowledgebaseEvent, T : SomaticEvent {
            return KnowledgebaseEventReader(source, readers.toList())
        }
    }

    override fun read(event: R): List<T> {
        val events = readers.map { Pair(it.javaClass.simpleName, it.read(event)) }.filterNot { it.second.isEmpty() }
        if (events.size > 1) {
            logger.warn("More than 1 reader (${events.joinToString(", ") { it.first }}) returned events for record: $event")
        }
        if (events.isEmpty()) {
            logger.warn("Could not extract somatic event from:\t$source\t${event.gene}\t${event.variant}")
        }
        return events.map { it.second }.flatten()
    }
}
