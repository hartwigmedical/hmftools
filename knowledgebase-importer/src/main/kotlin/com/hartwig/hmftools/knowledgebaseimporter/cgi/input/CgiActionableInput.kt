package com.hartwig.hmftools.knowledgebaseimporter.cgi.input

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CorrectedInput

data class CgiActionableInput(@get:JvmName("getGene_") private val Gene: String, override val transcript: String?, val Alteration: String,
                              val `Alteration type`: String, val gDNA: String?, val cDNA: String?, val individual_mutation: String?,
                              val `Primary Tumor type`: String, val `Evidence level`: String, val Association: String, val Drug: String?,
                              val `Drug family`: String?) : CsvData, CorrectedInput<CgiActionableInput>, CgiInput {

    override val gene: String = Gene
    override val variant: String = Alteration.substringAfter(":")
    override val gdna: String = gDNA.orEmpty()

    override fun correct(): CgiActionableInput? {
        return when {
            Alteration.contains("RET__TPCN1")     -> null
            Association == "No Responsive"        -> null
            Alteration.contains("ABL1__BCR")      -> copy(Alteration = Alteration.replace("ABL1__BCR", "BCR__ABL1"))
            Alteration.contains("BRD4__C15orf55") -> copy(Gene = "NUTM1", Alteration = Alteration.replace("BRD4__C15orf55", "BRD4__NUTM1"))
            Alteration.contains("PDGFRA__FIP1L1") -> copy(Alteration = Alteration.replace("PDGFRA__FIP1L1", "FIP1L1__PDGFRA"))
            Alteration.contains("PDGFB__COL1A1")  -> copy(Alteration = Alteration.replace("PDGFB__COL1A1", "COL1A1__PDGFB"))
            else                                  -> this
        }
    }
}
