package com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class TransvarProteinMatcherTest : StringSpec() {

    init {
        "matches substitutions"{
            TransvarProteinMatcher.matches("p.Trp24Cys") shouldBe true
            TransvarProteinMatcher.matches("Trp24Cys") shouldBe true
            TransvarProteinMatcher.matches("p.Trp24Ter") shouldBe true
            TransvarProteinMatcher.matches("Trp24Ter") shouldBe true
            TransvarProteinMatcher.matches("p.Trp24*") shouldBe true
            TransvarProteinMatcher.matches("Trp24*") shouldBe true
            TransvarProteinMatcher.matches("p.*757L") shouldBe true
            TransvarProteinMatcher.matches("*757L") shouldBe true
        }

        "matches substitutions with mixed AA lengths"{
            TransvarProteinMatcher.matches("p.V600Glu") shouldBe true
            TransvarProteinMatcher.matches("V600Glu") shouldBe true
            TransvarProteinMatcher.matches("p.Val600E") shouldBe true
            TransvarProteinMatcher.matches("Val600E") shouldBe true
        }

        "matches deletions"{
            TransvarProteinMatcher.matches("p.Gln8del") shouldBe true
            TransvarProteinMatcher.matches("Gln8del") shouldBe true
            TransvarProteinMatcher.matches("p.*4del") shouldBe true
            TransvarProteinMatcher.matches("*4del") shouldBe true
            TransvarProteinMatcher.matches("p.Cys28_Met30del") shouldBe true
            TransvarProteinMatcher.matches("Cys28_Met30del") shouldBe true
            TransvarProteinMatcher.matches("534_536del") shouldBe true
            TransvarProteinMatcher.matches("533_534del") shouldBe true
            TransvarProteinMatcher.matches("p.K601delK") shouldBe true
            TransvarProteinMatcher.matches("K601delK") shouldBe true
            TransvarProteinMatcher.matches("p.D898_E901delDVYE") shouldBe true
            TransvarProteinMatcher.matches("D898_E901delDVYE") shouldBe true

        }

        "matches duplications"{
            TransvarProteinMatcher.matches("p.G10dup") shouldBe true
            TransvarProteinMatcher.matches("G10dup") shouldBe true
            TransvarProteinMatcher.matches("p.Gly4_Gln6dup") shouldBe true
            TransvarProteinMatcher.matches("Gly4_Gln6dup") shouldBe true
            TransvarProteinMatcher.matches("p.G10dupG") shouldBe true
            TransvarProteinMatcher.matches("G10dupG") shouldBe true
            TransvarProteinMatcher.matches("p.V777_S779dupVGS") shouldBe true
            TransvarProteinMatcher.matches("V777_S779dupVGS") shouldBe true
        }

        "matches insertions"{
            TransvarProteinMatcher.matches("p.His4_Gln5insAla") shouldBe true
            TransvarProteinMatcher.matches("His4_Gln5insAla") shouldBe true
            TransvarProteinMatcher.matches("p.Lys2_Met3insGlnSerLys") shouldBe true
            TransvarProteinMatcher.matches("Lys2_Met3insGlnSerLys") shouldBe true
            TransvarProteinMatcher.matches("p.W603_E604insDREYEYDLKW") shouldBe true
            TransvarProteinMatcher.matches("W603_E604insDREYEYDLKW") shouldBe true
        }

        "matches deletion-insertions"{
            TransvarProteinMatcher.matches("p.Cys28delinsTrp") shouldBe true
            TransvarProteinMatcher.matches("Cys28delinsTrp") shouldBe true
            TransvarProteinMatcher.matches("p.Cys28delinsTrpVal") shouldBe true
            TransvarProteinMatcher.matches("Cys28delinsTrpVal") shouldBe true
            TransvarProteinMatcher.matches("p.Cys28_Lys29delinsTrp") shouldBe true
            TransvarProteinMatcher.matches("Cys28_Lys29delinsTrp") shouldBe true
            TransvarProteinMatcher.matches("p.C420_A423delinsW") shouldBe true
            TransvarProteinMatcher.matches("C420_A423delinsW") shouldBe true
            TransvarProteinMatcher.matches("p.Pro578_Lys579delinsLeuTer") shouldBe true
            TransvarProteinMatcher.matches("Pro578_Lys579delinsLeuTer") shouldBe true
        }

        "matches frameshifts" {
            TransvarProteinMatcher.matches("p.Glu5fs") shouldBe true
            TransvarProteinMatcher.matches("Glu5fs") shouldBe true
            TransvarProteinMatcher.matches("p.Glu5ValfsTer5") shouldBe true
            TransvarProteinMatcher.matches("Glu5ValfsTer5") shouldBe true
            TransvarProteinMatcher.matches("p.Glu5Valfs*5") shouldBe true
            TransvarProteinMatcher.matches("Glu5Valfs*5") shouldBe true
            TransvarProteinMatcher.matches("p.Ile327Argfs*?") shouldBe true
            TransvarProteinMatcher.matches("Ile327Argfs*?") shouldBe true
            TransvarProteinMatcher.matches("R259fs*15") shouldBe true
        }

        "does not match codon mutations" {
            TransvarProteinMatcher.matches("V600") shouldBe false
            TransvarProteinMatcher.matches("Val600") shouldBe false
        }

        "does not match complex variants"{
            TransvarProteinMatcher.matches("Fusion + V600E") shouldBe false
            TransvarProteinMatcher.matches("Fusion + Val600E") shouldBe false
            TransvarProteinMatcher.matches("Fusion + V600") shouldBe false
            TransvarProteinMatcher.matches("Fusion + Val600") shouldBe false
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

        "does not match gene fusion"{
            TransvarProteinMatcher.contains("RANBP2-ALK") shouldBe false
        }

        "does not match invalid frameshift"{
            TransvarProteinMatcher.matches("K442Nfs*") shouldBe false
            TransvarProteinMatcher.matches("N1333Gfs*") shouldBe false
        }

        "does not match insertions without range"{
            TransvarProteinMatcher.matches("T574insTQLPYD") shouldBe false
        }

        "finds protein mutation in complex variant"{
            TransvarProteinMatcher.contains("Fusion + V600E") shouldBe true
            TransvarProteinMatcher.contains("Fusion + Val600Glu") shouldBe true
            TransvarProteinMatcher.contains("Y175fs (c.526delA)") shouldBe true
            TransvarProteinMatcher.contains("F76del (c.226_228delTTC)") shouldBe true
        }
    }
}
