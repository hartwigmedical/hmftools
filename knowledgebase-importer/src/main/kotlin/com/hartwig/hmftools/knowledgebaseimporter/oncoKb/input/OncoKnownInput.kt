package com.hartwig.hmftools.knowledgebaseimporter.oncoKb.input

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CorrectedInput

data class OncoKnownInput(private val Isoform: String?, @get:JvmName("getGene_") private val Gene: String, val Alteration: String,
                          val `Mutation Effect`: String, val Oncogenicity: String) : CsvData, CorrectedInput<OncoKnownInput>,
        OncoKbInput {

    val reference = "$Gene $Alteration"
    override val transcript = Isoform
    override val gene = Gene
    override val variant = Alteration

    override fun correct(): OncoKnownInput? {
        return when {
            Alteration.contains("IGH-NKX2") && Gene == "NKX2-1" -> copy(Alteration = Alteration.replace("IGH-NKX2", "IGH-NKX2-1"))
            Alteration.contains("ROS1-CD74")                    -> copy(Alteration = Alteration.replace("ROS1-CD74", "CD74-ROS1"))
            Alteration.contains("RET-CCDC6")                    -> copy(Alteration = Alteration.replace("RET-CCDC6", "CCDC6-RET"))
            Alteration.contains("FGFR1OP1-FGFR1")               -> copy(Alteration = Alteration.replace("FGFR1OP-FGFR1", "FGFR1OP1-FGFR1"))
            Alteration.contains("EP300-MOZ")                    -> copy(Alteration = Alteration.replace("EP300-MOZ", "KAT6A-EP300"))
            Alteration.contains("EP300-MLL")                    -> copy(Alteration = Alteration.replace("EP300-MLL", "KMT2A-EP300"))
            Alteration.contains("BRD4-NUT")                     -> copy(Alteration = Alteration.replace("BRD4-NUT", "BRD4-NUTM1"))
            Alteration.contains("CEP110-FGFR1")                 -> copy(Alteration = Alteration.replace("CEP110-FGFR1", "CNTRL-FGFR1"))
            Alteration.contains("FGFR2-KIAA1967")               -> copy(Alteration = Alteration.replace("FGFR2-KIAA1967", "FGFR2-CCAR2"))
            Alteration.contains("FIG-ROS1")                     -> copy(Alteration = Alteration.replace("FIG-ROS1", "GOPC-ROS1"))
            Alteration.contains("GPIAP1-PDGFRB")                -> copy(Alteration = Alteration.replace("GPIAP1-PDGFRB", "CAPRIN1-PDGFRB"))
            Alteration.contains("IGL-MYC")                      -> copy(Alteration = Alteration.replace("IGL-MYC", "IGLC6-MYC"))
            Alteration.contains("KIAA1509-PDGFRB")              -> copy(Alteration = Alteration.replace("KIAA1509-PDGFRB", "CCDC88C-PDGFRB"))
            Alteration.contains("MLL-TET1")                     -> copy(Alteration = Alteration.replace("MLL-TET1", "KMT2A-TET1"))
            Alteration.contains("NPM-ALK")                      -> copy(Alteration = Alteration.replace("NPM-ALK", "NPM1-ALK"))
            Alteration.contains("PAX8-PPAR")                    -> copy(Alteration = Alteration.replace("PAX8-PPAR", "PAX8-PPARA"))
            Alteration.contains("SEC16A1-NOTCH1")               -> copy(Alteration = Alteration.replace("SEC16A1-NOTCH1", "SEC16A-NOTCH1"))
            Alteration.contains("TEL-JAK2")                     -> copy(Alteration = Alteration.replace("TEL-JAK2", "ETV6-JAK2"))
            Alteration.contains("TRA-NKX2-1")                   -> copy(Alteration = Alteration.replace("TRA-NKX2-1", "TRAC-NKX2-1"))
            Alteration == "p61BRAF-V600E"                       -> copy(Alteration = "V600E/V600K")
            Alteration == "Delta-NTRK1"                         -> null
            else                                                -> this
        }
    }
}
