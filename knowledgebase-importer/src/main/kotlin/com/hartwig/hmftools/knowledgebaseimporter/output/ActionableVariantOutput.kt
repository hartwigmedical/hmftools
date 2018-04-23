package com.hartwig.hmftools.knowledgebaseimporter.output

import com.hartwig.hmftools.common.variant.SomaticVariant

data class ActionableVariantOutput(private val gene: String, private val variant: SomaticVariant,
                                   private val actionability: Actionability) {
    companion object {
        val header = listOf("gene", "chromosome", "position", "ref", "alt", "source", "drug", "cancerType", "level", "type", "significance")
    }

    val record: List<String> = listOf(gene,
                                      variant.chromosome(),
                                      variant.position().toString(),
                                      variant.ref(),
                                      variant.alt(),
                                      actionability.source,
                                      actionability.drug,
                                      actionability.cancerType,
                                      actionability.level,
                                      actionability.type,
                                      actionability.significance)
}
