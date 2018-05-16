package com.hartwig.hmftools.knowledgebaseimporter

import com.hartwig.hmftools.knowledgebaseimporter.output.FusionPair
import com.hartwig.hmftools.knowledgebaseimporter.output.PromiscuousGene
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

        "finds potential five gene" {
            isFiveGene("KMT2A", "KMT2A-MLLT3", "-") shouldBe true
        }

        "finds potential three gene"{
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

        "extracts fiveGene with approx threeGene target from complex variant"{
            fiveGene("ABL1", "BCR-ABL T334I", "-") shouldBe "BCR"
        }

        "extracts threeGene with approx threeGene target from complex variant"{
            threeGene("ABL1", "BCR-ABL T334I", "-") shouldBe "ABL1"
        }

        "extracts fiveGene with threeGene target from complex variant"{
            fiveGene("ALK", "EML4-ALK L1196M", "-") shouldBe "EML4"
        }

        "extracts threeGene with threeGene target from complex variant"{
            threeGene("ALK", "EML4-ALK L1196M", "-") shouldBe "ALK"
        }

        "extracts fusion from simple fusion variant"{
            extractFusion("NTRK1", "LMNA-NTRK1", "-") shouldBe FusionPair("LMNA", "NTRK1")
        }

        "extracts fusion from truncating fusion variant"{
            extractFusion("FOS", "TRUNCATING FUSION", "-") shouldBe PromiscuousGene("FOS")
        }

        "extracts fusion from rearrangement variant"{
            extractFusion("ROS1", "REARRANGEMENT", "-") shouldBe PromiscuousGene("ROS1")
        }

        "extracts fusion from complex fusion variant"{
            extractFusion("ALK", "ALK FUSION G1202R", "-") shouldBe PromiscuousGene("ALK")
        }

        "extracts fusion from fusion amplification variant"{
            extractFusion("ALK", "EML4-ALK AMPLIFICATION", "-") shouldBe FusionPair("EML4", "ALK")
        }

        "extracts fusion from promiscuous fusion variant"{
            extractFusion("DUX4", "DUX4 FUSIONS", "-") shouldBe PromiscuousGene("DUX4")
        }

        "extracts fusion from ? separated variant"{
            extractFusion("FGFR1", "ERLIN2?FGFR1 Fusion", "?") shouldBe FusionPair("ERLIN2", "FGFR1")
        }

        "extracts fusion from ` - ` separated variant"{
            extractFusion("FGFR3", "FGFR3 - BAIAP2L1 Fusion", " - ") shouldBe FusionPair("FGFR3", "BAIAP2L1")
        }

        "extracts fusion with similar gene names"{
            extractFusion("FGFR1", "FGFR10-FGFR1 Fusion", "-") shouldBe FusionPair("FGFR10", "FGFR1")
        }

        "extracts fusion with - in gene name"{
            extractFusion("NKX2-1", "IGH-NKX2-1 Fusion", "-") shouldBe FusionPair("IGH", "NKX2-1")
        }
    }
}
