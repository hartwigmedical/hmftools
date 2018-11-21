package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.CnvEvent
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class IclusionCnvReaderTest : StringSpec() {
    companion object {
        private val event = IclusionEvent("ERBB2", "ANY", "ENST0000000")
    }

    init {
        "can read amplification" {
            IclusionCnvReader.read(event.copy(variant = "AMPLIFICATION")) shouldBe listOf(CnvEvent.amplification(event.gene))
            IclusionCnvReader.read(event.copy(variant = "OVEREXPRESSION")) shouldBe listOf(CnvEvent.amplification(event.gene))
        }

        "can read deletion" {
        //    IclusionCnvReader.read(actionable.copy(variant = "Deletion")) shouldBe listOf(CnvEvent.deletion(actionable.gene))
        }

        "will only read exact matches" {
            IclusionCnvReader.read(event.copy(variant = "amplification")) shouldBe emptyList<CnvEvent>()
            IclusionCnvReader.read(event.copy(variant = "DEL")) shouldBe emptyList<CnvEvent>()
            IclusionCnvReader.read(event.copy(variant = "V600")) shouldBe emptyList<CnvEvent>()
        }
    }

}