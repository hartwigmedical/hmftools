package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.FusionReader
import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class IclusionFusionReaderTest : StringSpec() {
    companion object {
        private const val FUSION_SEPARATOR = "-"
        private val fusionReader = FusionReader(separators = setOf(FUSION_SEPARATOR))

        private val event = IclusionEvent("BRAF-PIK3 Fusion", "ANY", "ENST0000000")
    }

    init {
        "can read fusion" {
            IclusionFusionReader.read(event.copy(variant = "rearrangement")) shouldBe listOf(fusionReader.read(event.gene, event.variant))
            IclusionFusionReader.read(event.copy(variant = "fusions")) shouldBe listOf(fusionReader.read(event.gene, event.variant))
            IclusionFusionReader.read(event.copy(variant = "fusion")) shouldBe listOf(fusionReader.read(event.gene, event.variant))
        }
    }
}
