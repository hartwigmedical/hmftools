package com.hartwig.hmftools.knowledgebaseimporter.knowledgebases

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class ProteinModelTest : StringSpec() {

    init {
        "matches protein mutation" {
            ProteinModel.isProteinAlteration("V600E") shouldBe true
            ProteinModel.isProteinAlteration("Val600Glu") shouldBe true
        }

        "codon mutation does not match protein mutation"{
            ProteinModel.isProteinAlteration("V600") shouldBe false
            ProteinModel.isProteinAlteration("Val600") shouldBe false
        }

        "complex variant does not match protein mutation"{
            ProteinModel.isProteinAlteration("Fusion + V600E") shouldBe false
            ProteinModel.isProteinAlteration("Fusion + Val600E") shouldBe false
        }

        "matches codon mutation" {
            ProteinModel.isCodonAlteration("V600") shouldBe true
            ProteinModel.isCodonAlteration("Val600") shouldBe true
        }

        "protein mutation does not match codon mutation" {
            ProteinModel.isCodonAlteration("V600E") shouldBe false
            ProteinModel.isCodonAlteration("Val600Glu") shouldBe false
        }

        "complex variant does not match codon mutation"{
            ProteinModel.isCodonAlteration("Fusion + V600") shouldBe false
            ProteinModel.isCodonAlteration("Fusion + Val600") shouldBe false
        }

        "finds protein mutation in complex variant"{
            ProteinModel.hasProteinAlteration("Fusion + V600E") shouldBe true
            ProteinModel.hasProteinAlteration("Fusion + Val600Glu") shouldBe true
        }

        "finds codon mutation in complex variant"{
            ProteinModel.hasCodonAlteration("Fusion + V600") shouldBe true
            ProteinModel.hasCodonAlteration("Fusion + Val600") shouldBe true
        }

        "finds codon mutation in protein variant in non-strict mode"{
            ProteinModel.hasCodonAlteration("V600E", false) shouldBe true
            ProteinModel.hasCodonAlteration("Val600Glu", false) shouldBe true
        }

        "does not finds codon mutation in protein variant in strict mode"{
            ProteinModel.hasCodonAlteration("V600E") shouldBe false
            ProteinModel.hasCodonAlteration("Val600Glu") shouldBe false
        }

        "does not match invalid codon number"{
            ProteinModel.isProteinAlteration("V0E") shouldBe false
            ProteinModel.hasProteinAlteration("V0E") shouldBe false
            ProteinModel.isCodonAlteration("V0") shouldBe false
            ProteinModel.hasCodonAlteration("V0") shouldBe false
        }
    }
}
