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
        variant.contains(Regex("MLL-MLLT3")) && gene == "KMT2A" -> copy(variant = variant.replace("MLL-MLLT3", "KMT2A-MLLT3"))
        variant == "BRAF-CUL1"                                  -> null
        variant == "ZNF198-FGFR1"                               -> copy(variant = "ZMYM2-FGFR1")
        transcript == "ENST0000023170.2"                        -> copy(representative_transcript = "ENST00000231790.2")
        transcript == "ENST000002564742"                        -> copy(representative_transcript = "ENST00000256474.2")
        else                                                    -> this
    }
}
