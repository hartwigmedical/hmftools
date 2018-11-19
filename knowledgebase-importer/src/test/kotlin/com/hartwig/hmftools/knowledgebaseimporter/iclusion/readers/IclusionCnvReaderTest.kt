package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.output.CnvEvent
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class IclusionCnvReaderTest : StringSpec() {
    companion object {
        private val actionable = IclusionEvent("BRAF", "V600", "ENST0000000")
    }

    init {
        "can read amplification" {
         //   IclusionCnvReader.read(actionable.copy() ) shouldBe listOf(CnvEvent.amplification(actionable.gene))
        }

        "can read deletion" {
        //    IclusionCnvReader.read(actionable.copy(variant = "Deletion")) shouldBe listOf(CnvEvent.deletion(actionable.gene))
        }

        "will only read exact matches" {
//            IclusionCnvReader.read(actionable.copy(variant = "amplification")) shouldBe emptyList<CnvEvent>()
//            IclusionCnvReader.read(actionable.copy(variant = "amp")) shouldBe emptyList<CnvEvent>()
//            IclusionCnvReader.read(actionable.copy(variant = "overexpression")) shouldBe emptyList<CnvEvent>()
//            IclusionCnvReader.read(actionable.copy(variant = "deletion")) shouldBe emptyList<CnvEvent>()
//            IclusionCnvReader.read(actionable.copy(variant = "del")) shouldBe emptyList<CnvEvent>()
//            IclusionCnvReader.read(actionable.copy(variant = "loss-of-function")) shouldBe emptyList<CnvEvent>()

        }
    }

}