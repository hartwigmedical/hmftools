package com.hartwig.hmftools.knowledgebaseimporter.output

import com.hartwig.hmftools.common.variant.SomaticVariant
import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent

sealed class ActionableEvent : SomaticEvent {
    abstract fun eventString(): String
}

sealed class FusionEvent : ActionableEvent()

data class FusionPair(val fiveGene: String, val threeGene: String) : FusionEvent(), CsvData {
    override fun eventString(): String {
        return "$fiveGene - $threeGene fusion"
    }
}

data class PromiscuousGene(val gene: String) : FusionEvent(), CsvData {
    override fun eventString(): String {
        return "$gene fusions"
    }
}

data class CnvEvent(val gene: String, val cnvType: String) : ActionableEvent(), CsvData {
    override fun eventString(): String {
        return "$gene $cnvType"
    }
}

data class SomaticVariantEvent(val gene: String, val chromosome: String, val position: String, val ref: String, val alt: String) :
        ActionableEvent(), CsvData {
    companion object {
        operator fun invoke(gene: String, variant: SomaticVariant): SomaticVariantEvent {
            return SomaticVariantEvent(gene, variant.chromosome(), variant.position().toString(), variant.ref(), variant.alt())
        }
    }

    override fun eventString(): String {
        return "$gene $chromosome:$position $ref->$alt"
    }
}

data class GenomicRangeEvent(val gene: String, val mutationTranscript: String, val chromosome: String, val start: String,
                             val stop: String, val geneTranscript: String) :
        ActionableEvent(), CsvData {

    override fun eventString(): String {
        return "$gene($mutationTranscript) $chromosome:$start-$stop"
    }
}
