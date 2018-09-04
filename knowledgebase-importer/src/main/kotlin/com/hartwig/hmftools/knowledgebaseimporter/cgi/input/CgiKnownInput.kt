package com.hartwig.hmftools.knowledgebaseimporter.cgi.input

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CorrectedInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.KnowledgebaseEvent

data class CgiKnownInput(override val gene: String, override val transcript: String, val context: String, val gdna: String,
                         val protein: String) : CsvData, CorrectedInput<CgiKnownInput>, KnowledgebaseEvent {
    override val variant: String = protein

    override fun correct(): CgiKnownInput? {
        return if (context != "somatic") null
        else this
    }
}
