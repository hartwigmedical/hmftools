package com.hartwig.hmftools.knowledgebaseimporter.civic.input

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CorrectedInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.KnowledgebaseEvent

data class CivicVariantInput(override val gene: String, private val representative_transcript: String, val variant_id: String,
                             override val variant: String, val chromosome: String, val start: String, val stop: String,
                             val reference_bases: String, val variant_bases: String, val hgvs_expressions: String,
                             private val variant_types: String) : CsvData,
        CorrectedInput<CivicVariantInput>, KnowledgebaseEvent {
    companion object {
        private const val RANGE_VARIANTS = "gene_variant|transcript_variant|exon_variant|coding_sequence_variant|protein_altering_variant"
    }

    override val transcript = representative_transcript
    val variantTypes: List<String> =
            if (variant_types.matches("N/A".toRegex(RegexOption.IGNORE_CASE))) emptyList()
            else variant_types.split(",").filterNot { it.isBlank() }
    val hasPosition = chromosome.isNotBlank() && start.isNotBlank() && stop.isNotBlank()
    val hasRefOrAlt = reference_bases.isNotBlank() || variant_bases.isNotBlank()
    private val isGenericMutation = variant.toLowerCase() == "mutation" || variantTypes.any { it.contains(RANGE_VARIANTS.toRegex()) }
    private val isGenericMissense = !variant.contains("+") && variantTypes.size == 1 && !variant.toLowerCase().contains(" and ") &&
            variantTypes.first() == "missense_variant"
    val isFusion = variantTypes.isNotEmpty() && variantTypes.all { it.contains("fusion") }
    val isCNV = variant == "AMPLIFICATION" || variant == "DELETION"
    val hgvs: String?
        get() {
            val matchResult = Regex("(ENST[0-9]+\\.[0-9+]:c\\.[0-9][^,\\t\\s\\n]+)").find(hgvs_expressions)
            return matchResult?.groupValues?.get(1)
        }
    val hasHgvs = hgvs != null
    val hasKnownVariant = hasPosition && hasRefOrAlt
    val hasVariant = hasKnownVariant || hasHgvs
    val isVariantRecord = variantTypes.none { it.contains("fusion") } && hasVariant
    val isGenomicRangeMutation = hasPosition && !hasRefOrAlt && (isGenericMutation || isGenericMissense)

    override fun correct(): CivicVariantInput? = when {
        variant.contains(Regex("MLL-MLLT3")) && gene == "KMT2A" -> copy(variant = variant.replace("MLL-MLLT3", "KMT2A-MLLT3"))
        variant == "BRAF-CUL1"                                  -> null
        else                                                    -> this
    }
}
