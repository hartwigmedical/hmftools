package com.hartwig.hmftools.knowledgebaseimporter.oncoKb.readers

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.KnowledgebaseAnyEventReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.KnowledgebaseEventReader
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoKnownInput

data class OncoAnyKnownEventReader<out T : SomaticEvent>(
        override val nestedEventReader: KnowledgebaseEventReader<OncoKnownInput, T>) :
        KnowledgebaseAnyEventReader<OncoKnownInput, T> {

    override fun mapper(input: OncoKnownInput): List<OncoKnownInput> {
        return if (input.variant.contains("/")) input.variant.split("/").map { input.copy(Alteration = it.trim()) }
        else emptyList()
    }
}
