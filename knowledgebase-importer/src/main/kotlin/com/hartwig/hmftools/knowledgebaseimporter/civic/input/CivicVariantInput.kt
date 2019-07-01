package com.hartwig.hmftools.knowledgebaseimporter.civic.input

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CorrectedInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.KnowledgebaseEvent

data class CivicVariantInput(override val gene: String, private val representative_transcript: String, val variant_id: String,
                             override val variant: String, val chromosome: String, val start: String, val stop: String,
                             val reference_bases: String, val variant_bases: String, val hgvs_expressions: String,
                             private val variant_types: String) : CsvData,
        CorrectedInput<CivicVariantInput>, KnowledgebaseEvent {

    override val transcript = representative_transcript
    val variantTypes: List<String> =
            if (variant_types.matches("N/A".toRegex(RegexOption.IGNORE_CASE))) emptyList()
            else variant_types.split(",").filterNot { it.isBlank() }
    val hasPosition = chromosome.isNotBlank() && start.isNotBlank() && stop.isNotBlank()
    val hasRefOrAlt = reference_bases.isNotBlank() || variant_bases.isNotBlank()

    override fun correct(): CivicVariantInput? = when {
        variant == "BRAF-CUL1"                                  -> null
        variant.contains(Regex("MLL-MLLT3")) && gene == "KMT2A" -> copy(variant = variant.replace("MLL-MLLT3", "KMT2A-MLLT3"))
        variant.contains("ZNF198-FGFR1")                        -> copy(variant = variant.replace("ZNF198-FGFR1", "ZMYM2-FGFR1"))
        variant.contains("NTRK2-STRN")                          -> copy(variant = variant.replace("NTRK2-STRN", "STRN-NTRK2"))
        variant.contains("NPM-ALK")                             -> copy(variant = variant.replace("NPM-ALK", "NPM1-ALK"))
        variant.contains("WASFL-BRAF")                             -> copy(variant = variant.replace("WASFL-BRAF", "WASF1-BRAF"))
        representative_transcript == "ENST0000023170.2"         -> copy(representative_transcript = "ENST00000231790.2")
        representative_transcript == "ENST000002564742"         -> copy(representative_transcript = "ENST00000256474.2")
        else                                                    -> this
    }
}
