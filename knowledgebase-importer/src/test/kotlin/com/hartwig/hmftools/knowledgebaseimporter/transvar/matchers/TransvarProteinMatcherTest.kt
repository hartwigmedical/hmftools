package com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class TransvarProteinMatcherTest : StringSpec() {

    init {
        "matches protein mutation" {
            TransvarProteinMatcher.matches("V600E") shouldBe true
            TransvarProteinMatcher.matches("Val600Glu") shouldBe true
        }

        "complex variant does not match"{
            TransvarProteinMatcher.matches("Fusion + V600E") shouldBe false
            TransvarProteinMatcher.matches("Fusion + Val600E") shouldBe false
            TransvarProteinMatcher.matches("Fusion + V600") shouldBe false
            TransvarProteinMatcher.matches("Fusion + Val600") shouldBe false
        }

        "matches codon mutation" {
            TransvarProteinMatcher.matches("V600") shouldBe true
            TransvarProteinMatcher.matches("Val600") shouldBe true
        }

        "finds protein mutation in complex variant"{
            TransvarProteinMatcher.contains("Fusion + V600E") shouldBe true
            TransvarProteinMatcher.contains("Fusion + Val600Glu") shouldBe true
            TransvarProteinMatcher.contains("Fusion + V600") shouldBe true
            TransvarProteinMatcher.contains("Fusion + Val600") shouldBe true
            TransvarProteinMatcher.contains("Y175fs (c.526delA)") shouldBe true
            TransvarProteinMatcher.contains("F76del (c.226_228delTTC)") shouldBe true
            TransvarProteinMatcher.contains("R108ins (c.324InsCGC)") shouldBe true
            TransvarProteinMatcher.contains("G12/G13") shouldBe true
        }

        "does not match other variants" {
            TransvarProteinMatcher.contains("DUX4 Fusions") shouldBe false
            TransvarProteinMatcher.contains("SH2 DOMAIN MUTATION") shouldBe false
            TransvarProteinMatcher.contains("DEL 485-490") shouldBe false
            TransvarProteinMatcher.contains("RARE EX 18-21 MUT") shouldBe false
            TransvarProteinMatcher.contains("Exon 12 splice site insertion") shouldBe false
            TransvarProteinMatcher.contains("EWSR1-FLI1 Type 1") shouldBe false
            TransvarProteinMatcher.contains("DPYD*2A HOMOZYGOSITY") shouldBe false
            TransvarProteinMatcher.contains("3` UTR Polymorphism") shouldBe false
            TransvarProteinMatcher.contains("Splicing alteration (c.463+2C>T)") shouldBe false
        }

        "codon mutation does not match gene fusion"{
            TransvarProteinMatcher.contains("RANBP2-ALK") shouldBe false
        }

        "matches hgvs protein mutation"{
            TransvarProteinMatcher.matches("Y175fs") shouldBe true
            TransvarProteinMatcher.matches("F76del") shouldBe true
            TransvarProteinMatcher.matches("R108ins") shouldBe true
            TransvarProteinMatcher.matches("Cys28delinsTrpVal") shouldBe true
            TransvarProteinMatcher.matches("Arg97GlyfsTer26") shouldBe true
            TransvarProteinMatcher.matches("Arg97Glyfs*26") shouldBe true
            TransvarProteinMatcher.matches("Arg97_Leu833delinsGlyfsTer26") shouldBe true
            TransvarProteinMatcher.matches("Glu5Valfs*5") shouldBe true
            TransvarProteinMatcher.matches("Tyr4*") shouldBe true
            TransvarProteinMatcher.matches("Asp2Metfs*4") shouldBe true
            TransvarProteinMatcher.matches("Glu5Valfs*5") shouldBe true
            TransvarProteinMatcher.matches("Ile327fs") shouldBe true
            TransvarProteinMatcher.matches("Gln151Thrfs*9") shouldBe true
        }

        "matches hgvs protein substitutions"{
            TransvarProteinMatcher.matches("*757L") shouldBe true
            TransvarProteinMatcher.matches("Trp26Cys") shouldBe true
            TransvarProteinMatcher.matches("Phe2_Met46del") shouldBe true
            TransvarProteinMatcher.matches("Met1?") shouldBe true
            TransvarProteinMatcher.matches("Trp26*") shouldBe true
            TransvarProteinMatcher.matches("Trp26_Leu833del") shouldBe true
        }

        "matches hgvs protein deletions"{
            TransvarProteinMatcher.matches("Gln8del") shouldBe true
            TransvarProteinMatcher.matches("Cys28_Met30del") shouldBe true
            TransvarProteinMatcher.matches("Met1_Leu833del") shouldBe true
            TransvarProteinMatcher.matches("Met1_Lys45del") shouldBe true
            TransvarProteinMatcher.matches("W26*") shouldBe true
        }

        "matches hgvs protein duplications"{
            TransvarProteinMatcher.matches("G10dup") shouldBe true
            TransvarProteinMatcher.matches("T1151dup") shouldBe true
            TransvarProteinMatcher.matches("Gly4_Gln6dup") shouldBe true
            TransvarProteinMatcher.matches("His7_Gln8dup") shouldBe true
            TransvarProteinMatcher.matches("A502_Y503dup") shouldBe true
        }

        "matches hgvs protein insertions"{
            TransvarProteinMatcher.matches("Lys2_Met3insGlnSerLys") shouldBe true
            TransvarProteinMatcher.matches("Pro2_Ile3insGlyTer") shouldBe true
            TransvarProteinMatcher.matches("Ile3_Ile3418delinsGly") shouldBe true
            TransvarProteinMatcher.matches("W603_E604insDREYEYDLKW") shouldBe true
        }

        "matches hgvs protein indels"{
            TransvarProteinMatcher.matches("Cys28_Lys29delinsTrp") shouldBe true
            TransvarProteinMatcher.matches("Cys28delinsTrpVal") shouldBe true
            TransvarProteinMatcher.matches("Pro578_Lys579delinsLeuTer") shouldBe true
            TransvarProteinMatcher.matches("Pro578_Gln598de") shouldBe true
            TransvarProteinMatcher.matches("Arg97Glyfs*26") shouldBe true
            TransvarProteinMatcher.matches("Arg97GlyfsTer26") shouldBe true
            TransvarProteinMatcher.matches("Arg97fs") shouldBe true
            TransvarProteinMatcher.matches("Arg97_Leu833delinsGlyfsTer26") shouldBe true
            TransvarProteinMatcher.matches("p.C420_A423delinsW") shouldBe true
            TransvarProteinMatcher.matches("Arg97fs") shouldBe true
            TransvarProteinMatcher.matches("Arg97Gly") shouldBe true
            TransvarProteinMatcher.matches("Arg97ProfsTer23") shouldBe true
            TransvarProteinMatcher.matches("Arg97Profs*23") shouldBe true
            TransvarProteinMatcher.matches("Glu5Valfs*5") shouldBe true
            TransvarProteinMatcher.matches("Tyr4*") shouldBe true
            TransvarProteinMatcher.matches("Asp2Metfs*4") shouldBe true
            TransvarProteinMatcher.matches("Asp2fs") shouldBe true
            TransvarProteinMatcher.matches("Glu5Valfs*5") shouldBe true
            TransvarProteinMatcher.matches("Glu5fs") shouldBe true
            TransvarProteinMatcher.matches("Ile327Argfs*?") shouldBe true
            TransvarProteinMatcher.matches("Ile327fs") shouldBe true
            TransvarProteinMatcher.matches("Gln151Thrfs*9") shouldBe true
        }
    }
}
