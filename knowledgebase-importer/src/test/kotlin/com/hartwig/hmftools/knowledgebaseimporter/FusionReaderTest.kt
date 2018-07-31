package com.hartwig.hmftools.knowledgebaseimporter

import com.hartwig.hmftools.knowledgebaseimporter.output.FusionPair
import com.hartwig.hmftools.knowledgebaseimporter.output.PromiscuousGene
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class FusionReaderTest : StringSpec() {

    init {
        "finds five gene" {
            FusionReader.isFiveGene("KMT2A", "KMT2A-MLLT3", "-") shouldBe true
        }

        "finds three gene"{
            FusionReader.isThreeGene("KMT2A", "MLLT3-KMT2A", "-") shouldBe true
        }

        "finds potential five gene" {
            FusionReader.isFiveGene("KMT2A", "KMT2A-MLLT3", "-") shouldBe true
        }

        "finds potential three gene"{
            FusionReader.isThreeGene("KMT2A", "MLLT3-KMT2A", "-") shouldBe true
        }

        "finds five gene with dash"{
            FusionReader.isFiveGene("NKX2-1", "NKX2-1-IGH", "-") shouldBe true
        }

        "finds three gene with dash"{
            FusionReader.isThreeGene("NKX2-1", "IGH-NKX2-1", "-") shouldBe true
        }

        "extracts fiveGene with fiveGene target"{
            FusionReader.fiveGene("KMT2A", "KMT2A-MLLT3", "-") shouldBe "KMT2A"
        }

        "extracts threeGene with fiveGene target"{
            FusionReader.threeGene("KMT2A", "KMT2A-MLLT3", "-") shouldBe "MLLT3"
        }

        "extracts fiveGene with threeGene target"{
            FusionReader.fiveGene("MLLT3", "KMT2A-MLLT3", "-") shouldBe "KMT2A"
        }

        "extracts threeGene with threeGene target"{
            FusionReader.threeGene("MLLT3", "KMT2A-MLLT3", "-") shouldBe "MLLT3"
        }

        "extracts fiveGene with approx threeGene target from complex variant"{
            FusionReader.fiveGene("ABL1", "BCR-ABL T334I", "-") shouldBe "BCR"
        }

        "extracts threeGene with approx threeGene target from complex variant"{
            FusionReader.threeGene("ABL1", "BCR-ABL T334I", "-") shouldBe "ABL1"
        }

        "extracts fiveGene with threeGene target from complex variant"{
            FusionReader.fiveGene("ALK", "EML4-ALK L1196M", "-") shouldBe "EML4"
        }

        "extracts threeGene with threeGene target from complex variant"{
            FusionReader.threeGene("ALK", "EML4-ALK L1196M", "-") shouldBe "ALK"
        }

        "extracts fusion from simple fusion variant"{
            FusionReader.extractFusion("NTRK1", "LMNA-NTRK1", "-") shouldBe
                    FusionPair("LMNA", "NTRK1")
        }

        "extracts fusion from truncating fusion variant"{
            FusionReader.extractFusion("FOS", "TRUNCATING FUSION", "-") shouldBe PromiscuousGene("FOS")
        }

        "extracts fusion from rearrangement variant"{
            FusionReader.extractFusion("ROS1", "REARRANGEMENT", "-") shouldBe PromiscuousGene("ROS1")
        }

        "extracts fusion from complex fusion variant"{
            FusionReader.extractFusion("ALK", "ALK FUSION G1202R", "-") shouldBe PromiscuousGene("ALK")
        }

        "extracts fusion from fusion amplification variant"{
            FusionReader.extractFusion("ALK", "EML4-ALK AMPLIFICATION", "-") shouldBe
                    FusionPair("EML4", "ALK")
        }

        "extracts fusion from promiscuous fusion variant"{
            FusionReader.extractFusion("DUX4", "DUX4 FUSIONS", "-") shouldBe PromiscuousGene("DUX4")
        }

        "extracts fusion from ? separated variant"{
            FusionReader.extractFusion("FGFR1", "ERLIN2?FGFR1 Fusion", "?") shouldBe
                    FusionPair("ERLIN2", "FGFR1")
        }

        "extracts fusion from ` - ` separated variant"{
            FusionReader.extractFusion("FGFR3", "FGFR3 - BAIAP2L1 Fusion", " - ") shouldBe
                    FusionPair("FGFR3", "BAIAP2L1")
        }

        "extracts fusion with similar gene names"{
            FusionReader.extractFusion("FGFR1", "FGFR10-FGFR1 Fusion", "-") shouldBe
                    FusionPair("FGFR10", "FGFR1")
        }

        "extracts fusion with - in gene name"{
            FusionReader.extractFusion("NKX2-1", "IGH-NKX2-1 Fusion", "-") shouldBe
                    FusionPair("IGH", "NKX2-1")
        }
    }
}
