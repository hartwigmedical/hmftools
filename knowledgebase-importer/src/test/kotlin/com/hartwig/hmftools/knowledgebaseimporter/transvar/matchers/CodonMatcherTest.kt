package com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class CodonMatcherTest : StringSpec() {

    init {
        "matches codon variant"{
            CodonMatcher.matches("V600") shouldBe true
            CodonMatcher.matches("Val600") shouldBe true
        }

        "finds protein mutation in complex variant"{
            CodonMatcher.contains("Fusion + V600") shouldBe true
            CodonMatcher.contains("Fusion + Val600") shouldBe true
            CodonMatcher.contains("G12/G13") shouldBe true
        }
    }
}
