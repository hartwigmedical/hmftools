package com.hartwig.hmftools.knowledgebaseimporter.oncoKb.readers

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.KnowledgebaseAnyEventReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.KnowledgebaseEventReader
import com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input.OncoActionableInput

data class OncoAnyActionableEventReader<out T : SomaticEvent>(
        override val nestedEventReader: KnowledgebaseEventReader<OncoActionableInput, T>) :
        KnowledgebaseAnyEventReader<OncoActionableInput, T> {

    override fun mapper(input: OncoActionableInput): List<OncoActionableInput> {
        return if (input.variant.contains("/")) input.variant.split("/").map { input.copy(Alteration = it.trim()) }
        else emptyList()
    }
}
