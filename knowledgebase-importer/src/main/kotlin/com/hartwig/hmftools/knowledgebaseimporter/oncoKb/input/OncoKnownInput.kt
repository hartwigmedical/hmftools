package com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CorrectedInput

data class OncoKnownInput(private val Isoform: String?, val Gene: String, val Alteration: String, val `Mutation Effect`: String,
                          val Oncogenicity: String) : CsvData, CorrectedInput<OncoKnownInput> {

    val reference = "$Gene $Alteration"
    val transcript = Isoform.orEmpty()

    override fun correct(): OncoKnownInput {
        return when {
            Alteration.contains(Regex("IGH-NKX2")) && Gene == "NKX2-1" -> copy(Alteration = Alteration.replace("IGH-NKX2", "IGH-NKX2-1"))
            else                                                       -> this
        }
    }
}
