package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CodonRangeMutations
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class IclusionCodonRangeReaderTest : StringSpec() {
    companion object {
        private val event = IclusionEvent("NRAS", "G12-G13", "ENST0000000")
    }

    init {
        "can read G12-G13" {
            IclusionCodonRangeReader.read(event) shouldBe listOf(CodonRangeMutations(event.gene, event.transcript, 12, 13))
        }
    }
}
