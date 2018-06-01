package com.hartwig.hmftools.knowledgebaseimporter.output

import com.hartwig.hmftools.common.variant.SomaticVariant
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent

sealed class ActionableEvent : SomaticEvent

sealed class FusionEvent : ActionableEvent()

data class FusionPair(val fiveGene: String, val threeGene: String) : FusionEvent() {
    companion object {
        val header = listOf("fiveGene", "threeGene")
    }

    val record: List<String> = listOf(fiveGene, threeGene)

    override fun toString(): String {
        return "$fiveGene - $threeGene fusion"
    }
}

data class PromiscuousGene(val gene: String) : FusionEvent() {
    companion object {
        val header = listOf("gene")
    }

    val record: List<String> = listOf(gene)

    override fun toString(): String {
        return "$gene fusions"
    }
}

data class CnvEvent(val gene: String, val cnvType: String) : ActionableEvent() {
    companion object {
        val header = listOf("gene", "cnvType")
    }

    val record: List<String> = listOf(gene, cnvType)

    override fun toString(): String {
        return "$gene $cnvType"
    }
}

data class SomaticVariantEvent(val gene: String, val chromosome: String, val position: String, val ref: String, val alt: String) :
        ActionableEvent() {
    companion object {
        val header = listOf("chromosome", "position", "ref", "alt")

        operator fun invoke(gene: String, variant: SomaticVariant): SomaticVariantEvent {
            return SomaticVariantEvent(gene, variant.chromosome(), variant.position().toString(), variant.ref(), variant.alt())
        }
    }

    val record: List<String> = listOf(chromosome, position, ref, alt)

    override fun toString(): String {
        return "$gene $chromosome:$position $ref->$alt"
    }
}
