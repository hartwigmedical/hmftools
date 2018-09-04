package com.hartwig.hmftools.knowledgebaseimporter.cgi.input

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CorrectedInput

data class CgiKnownInput(override val gene: String, override val transcript: String, val context: String, override val gdna: String,
                         val protein: String) : CsvData, CorrectedInput<CgiKnownInput>, CgiInput {
    override val variant: String = protein

    override fun correct(): CgiKnownInput? {
        return if (context != "somatic") null
        else this
    }
}
