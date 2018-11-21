package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.GeneMutations
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class IclusionGeneMutationReaderTest : StringSpec() {
    companion object {
        private val event = IclusionEvent("NRAS", "ANY", "ENST0000000")
    }

    init {
        "can match gene mutation " {
            IclusionGeneMutationReader.read(event.copy(variant = "MUTATION")) shouldBe listOf(GeneMutations(event.gene, event.transcript))
            IclusionGeneMutationReader.read(event.copy(variant = "ALTERATION")) shouldBe listOf(GeneMutations(event.gene, event.transcript))
            IclusionGeneMutationReader.read(event.copy(variant = "ACTIVATING MUTATION")) shouldBe listOf(GeneMutations(event.gene, event.transcript))
            IclusionGeneMutationReader.read(event.copy(variant = "INACTIVATING MUTATION")) shouldBe listOf(GeneMutations(event.gene, event.transcript))
        }

        "can not match gene mutation " {
            IclusionGeneMutationReader.read(event.copy(variant = "MUT")) shouldBe emptyList<GeneMutations>()
        }
    }
}
