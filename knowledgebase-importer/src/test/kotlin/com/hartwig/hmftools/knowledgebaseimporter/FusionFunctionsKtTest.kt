package com.hartwig.hmftools.knowledgebaseimporter

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class FusionFunctionsKtTest : StringSpec() {
    init {
        "finds five gene" {
            isFiveGene("KMT2A", "KMT2A-MLLT3", "-") shouldBe true
        }
        "finds three gene"{
            isThreeGene("KMT2A", "MLLT3-KMT2A", "-") shouldBe true
        }
        "finds five gene with dash"{
            isFiveGene("NKX2-1", "NKX2-1-IGH", "-") shouldBe true
        }
        "finds three gene with dash"{
            isThreeGene("NKX2-1", "IGH-NKX2-1", "-") shouldBe true
        }
        "extracts fiveGene with fiveGene target"{
            fiveGene("KMT2A", "KMT2A-MLLT3", "-") shouldBe "KMT2A"
        }
        "extracts threeGene with fiveGene target"{
            threeGene("KMT2A", "KMT2A-MLLT3", "-") shouldBe "MLLT3"
        }
        "extracts fiveGene with threeGene target"{
            fiveGene("MLLT3", "KMT2A-MLLT3", "-") shouldBe "KMT2A"
        }
        "extracts threeGene with threeGene target"{
            threeGene("MLLT3", "KMT2A-MLLT3", "-") shouldBe "MLLT3"
        }
    }
}
