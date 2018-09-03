package com.hartwig.hmftools.knowledgebaseimporter.iclusion.readers

import com.hartwig.hmftools.knowledgebaseimporter.iclusion.IclusionEvent
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.CodonMutations
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class IclusionCodonReaderTest : StringSpec() {
    companion object {
        private val codonEvent = IclusionEvent("BRAF", "V600", "ENST0000000")
        private val expectedVariant = CodonMutations(codonEvent.gene, codonEvent.transcript, 600)
    }

    init {
        "matches codon variant"{
            IclusionCodonReader.read(codonEvent) shouldBe listOf(expectedVariant)
            IclusionCodonReader.read(codonEvent.copy(variant = "Val600")) shouldBe listOf(expectedVariant)
        }

        "does not match protein variant"{
            IclusionCodonReader.read(codonEvent.copy(variant = "V600E")) shouldBe emptyList<SomaticEvent>()
        }

        "does not match multi-AA variant"{
            IclusionCodonReader.read(codonEvent.copy(variant = "VV600E")) shouldBe emptyList<SomaticEvent>()
        }
    }
}
