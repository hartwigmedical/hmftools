package com.hartwig.hmftools.knowledgebaseimporter.cosmic.input

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CorrectedInput

class CosmicKnownInput{}

//data class CosmicKnownInput(@get:JvmName("getGene_") private val gene: String, val Alteration: String) :
//        CsvData, CorrectedInput<CosmicKnownInput>,
//        ComicInput {
//
//
//        override fun correct(): CosmicKnownInput? {
//            return when {
//                Alteration.equals("CIC-DUX4L1") -> copy(Alteration = Alteration.replace("CIC-DUX4L1", "CIC-DUX4"))
//                Alteration.equals("CTAGE5-SIP1") -> copy(Alteration = Alteration.replace("CTAGE5-SIP1", "CTAGE5-GEMIN2"))
//                Alteration.equals("KMT2A-KIAA0284") -> copy(Alteration = Alteration.replace("KMT2A-KIAA0284", "KMT2A-CEP170B"))
//                Alteration.equals("NF1-ACCN1") -> copy(Alteration = Alteration.replace("NF1-ACCN1", "NF1-ASIC2"))
//                Alteration.equals("YWHAE-FAM22A") -> copy(Alteration = Alteration.replace("YWHAE-FAM22A", "YWHAE-NUTM2A"))
//                Alteration.equals("NF1-ACCN1") -> copy(Alteration = Alteration.replace("NF1-ACCN1", "NF1-ASIC2"))
//
//                else -> this
//            }
//        }
//    }