package com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfLevel
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfResponse

data class OncoActionableInput(private val Isoform: String?, val Gene: String, val Alteration: String, val `Cancer Type`: String,
                               val `Drugs(s)`: String, private val Level: String) : CsvData {

    val level = if (Level.startsWith("R")) Level.drop(1) else Level
    val hmfLevel = HmfLevel(Level)
    val significance = if (Level.startsWith("R")) HmfResponse.Resistant else HmfResponse.Responsive
    val reference = "$Gene $Alteration"
    val transcript = Isoform.orEmpty()
}
