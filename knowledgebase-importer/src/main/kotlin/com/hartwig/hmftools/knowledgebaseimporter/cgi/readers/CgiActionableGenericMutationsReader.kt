package com.hartwig.hmftools.knowledgebaseimporter.cgi.readers

import com.hartwig.hmftools.knowledgebaseimporter.cgi.input.CgiActionableInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.GenericMutation
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.readers.SomaticEventReader

object CgiActionableGenericMutationsReader : SomaticEventReader<CgiActionableInput, GenericMutation> {
    private val readers = listOf(CgiActionableCodonMutationsReader, CgiActionableCodonRangeMutationsReader,
                                 CgiActionableGeneMutationsReader)

    override fun read(event: CgiActionableInput): List<GenericMutation> {
        if (event.`Alteration type` == "MUT") {
            return readAlterations(event).map { event.copy(Alteration = it) }.flatMap { alteration -> readers.flatMap { it.read(alteration) } }
        }
        return emptyList()
    }

    private fun readAlterations(input: CgiActionableInput) = input.Alteration.substringAfter(":").split(",").map { it.trim() }

}