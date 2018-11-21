package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.ExonMutations
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class IclusionExonMutationReaderTest : StringSpec() {
    companion object {
        private val event = IclusionEvent("KRAS", "exon 100 mutation", "ENST0000000")
        private val eventRange = IclusionEvent("KRAS", "exon 1-3 mutation", "ENST0000000")

    }

    init {
        "can read exon mutation" {
            IclusionExonMutationReader.read(event) shouldBe listOf(ExonMutations(event.gene, event.transcript, 100))
        }

        "can read exon mutation range" {
            IclusionExonMutationReader.read(eventRange) shouldBe listOf(ExonMutations(eventRange.gene, eventRange.transcript, 1),
                    ExonMutations(eventRange.gene, eventRange.transcript, 2), ExonMutations(eventRange.gene, eventRange.transcript, 3))
        }
    }
}