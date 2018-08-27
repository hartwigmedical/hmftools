package com.hartwig.hmftools.knowledgebaseimporter.civic.input

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CorrectedInput

data class CivicVariantInput(val gene: String, val representative_transcript: String, val variant_id: String, val variant: String,
                             val chromosome: String, val start: String, val stop: String, val reference_bases: String,
                             val variant_bases: String, val hgvs_expressions: String, private val variant_types: String) : CsvData,
        CorrectedInput<CivicVariantInput> {
    companion object {
        private const val RANGE_VARIANTS = "gene_variant|transcript_variant|exon_variant|coding_sequence_variant|protein_altering_variant"
    }

    val variantTypes: List<String> =
            if (variant_types.matches("N/A".toRegex(RegexOption.IGNORE_CASE))) emptyList()
            else variant_types.split(",").filterNot { it.isBlank() }
    private val hasPosition = chromosome.isNotBlank() && start.isNotBlank() && stop.isNotBlank()
    private val hasRefOrAlt = reference_bases.isNotBlank() || variant_bases.isNotBlank()
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

    override fun correct(): CivicVariantInput = when {
        variant.contains(Regex("MLL-MLLT3")) && gene == "KMT2A" -> copy(variant = variant.replace("MLL-MLLT3", "KMT2A-MLLT3"))
        else                                                    -> this
    }
}
