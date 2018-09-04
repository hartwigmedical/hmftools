package com.hartwig.hmftools.knowledgebaseimporter.cgi.input

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CorrectedInput

data class CgiActionableInput(@get:JvmName("getGene_") private val Gene: String, override val transcript: String?, val Alteration: String,
                              val `Alteration type`: String, val gDNA: String?, val cDNA: String?, val individual_mutation: String?,
                              val `Primary Tumor type`: String, val `Evidence level`: String, val Association: String, val Drug: String?,
                              val `Drug family`: String?, override val variant: String = "") : CsvData, CorrectedInput<CgiActionableInput>,
        CgiInput {

    override val gene: String = Gene
    override val gdna: String = gDNA.orEmpty()

    override fun correct(): CgiActionableInput? {
        return when {
            Alteration.contains("RET-TPCN1")     -> null
            Association == "No Responsive"       -> null
            Alteration.contains("ABL1-BCR")      -> copy(Alteration = Alteration.replace("ABL1-BCR", "BCR-ABL1"))
            Alteration.contains("PDGFRA-FIP1L1") -> copy(Alteration = Alteration.replace("PDGFRA-FIP1L1", "FIP1L1-PDGFRA"))
            Alteration.contains("PDGFB-COL1A1")  -> copy(Alteration = Alteration.replace("PDGFB-COL1A1", "COL1A1-PDGFB"))
            else                                 -> this
        }
    }
}
