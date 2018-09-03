package com.hartwig.hmftools.knowledgebaseimporter.civic.readers

import com.hartwig.hmftools.knowledgebaseimporter.civic.input.CivicVariantInput
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.KnowledgebaseVariant
import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.SomaticEvent
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class CivicKnowledgebaseVariantReaderTest : StringSpec() {
    companion object {
        val variant = CivicVariantInput("BRAF", "ENST00000288602", "", "V600E", "7",
                                        "140453136", "140453136", "A", "T", "", "")
        private val expectedResult = KnowledgebaseVariant("BRAF", "7", 140453136, "A", "T")
    }

    init {
        "can read knowledgebase substitution variant" {
            CivicKnowledgebaseVariantReader.read(variant) shouldBe listOf(expectedResult)
        }

        "can read knowledgebase insertion variant" {
            CivicKnowledgebaseVariantReader.read(variant.copy(reference_bases = "")) shouldBe listOf(expectedResult.copy(ref = null))
        }

        "can read knowledgebase deletion variant" {
            CivicKnowledgebaseVariantReader.read(variant.copy(variant_bases = "")) shouldBe listOf(expectedResult.copy(alt = null))
        }

        "does not read variant when chromosome is missing" {
            CivicKnowledgebaseVariantReader.read(variant.copy(chromosome = "")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read variant when start is missing" {
            CivicKnowledgebaseVariantReader.read(variant.copy(start = "")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read variant when stop is missing" {
            CivicKnowledgebaseVariantReader.read(variant.copy(stop = "")) shouldBe emptyList<SomaticEvent>()
        }

        "does not read variant when both ref and alt are missing" {
            CivicKnowledgebaseVariantReader.read(variant.copy(reference_bases = "", variant_bases = "")) shouldBe emptyList<SomaticEvent>()
        }
    }
}
