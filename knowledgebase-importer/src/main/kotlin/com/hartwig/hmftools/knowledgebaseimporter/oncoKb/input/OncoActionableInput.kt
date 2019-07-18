package com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CorrectedInput
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfLevel
import com.hartwig.hmftools.knowledgebaseimporter.output.HmfResponse

data class OncoActionableInput(private val Isoform: String?, @get:JvmName("getGene_") private val Gene: String, val Alteration: String,
                               val `Cancer Type`: String, val `Drugs(s)`: String, private val Level: String) : CsvData,
        CorrectedInput<OncoActionableInput>, OncoKbInput {

    val level = if (Level.startsWith("R")) Level.drop(1) else Level
    val hmfLevel = HmfLevel(Level)
    val hmfResponse = if (Level.startsWith("R")) HmfResponse.Resistant else HmfResponse.Responsive
    val reference = Alteration
    override val transcript = Isoform.orEmpty()
    override val gene: String = Gene
    override val variant: String = Alteration

    override fun correct(): OncoActionableInput {
        val cancerType = if (`Cancer Type` == "Melanoma") "Skin Melanoma" else `Cancer Type`
        val alteration = if (Alteration == "p61BRAF-V600E") "V600E/V600K" else Alteration
        return copy(`Cancer Type` = cancerType, Alteration = alteration)
    }
}
