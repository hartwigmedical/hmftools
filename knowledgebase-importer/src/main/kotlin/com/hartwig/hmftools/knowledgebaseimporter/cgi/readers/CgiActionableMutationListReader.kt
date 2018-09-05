package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.KnowledgebaseAnyEventReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.KnowledgebaseEventReader
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.ProteinAnnotationReader

object CgiActionableMutationListReader : KnowledgebaseAnyEventReader<CgiActionableInput, SomaticEvent> {
    override val nestedEventReader = KnowledgebaseEventReader("cgi", ProteinAnnotationReader, CgiActionableCodonMutationsReader,
                                                              CgiActionableCodonRangeMutationsReader, CgiActionableGeneMutationsReader)

    override fun mapper(input: CgiActionableInput): List<CgiActionableInput> {
        return if (input.`Alteration type` == "MUT" && input.Alteration.contains(","))
            input.Alteration.substringAfter(":").split(",").map { input.copy(Alteration = it.trim()) }
        else emptyList()
    }
}
