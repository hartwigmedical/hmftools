package com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class TransvarCDnaMatcherTest : StringSpec() {

    init {
        "matches substitutions" {
            TransvarCDnaMatcher.matches("c.76A>C") shouldBe true
            TransvarCDnaMatcher.matches("76A>C") shouldBe true
            TransvarCDnaMatcher.matches("c.-14G>C") shouldBe true
            TransvarCDnaMatcher.matches("-14G>C") shouldBe true
            TransvarCDnaMatcher.matches("c.88+1G>T") shouldBe true
            TransvarCDnaMatcher.matches("88+1G>T") shouldBe true
            TransvarCDnaMatcher.matches("c.89-2A>C") shouldBe true
            TransvarCDnaMatcher.matches("89-2A>C") shouldBe true
            TransvarCDnaMatcher.matches("c.*46T>A") shouldBe true
            TransvarCDnaMatcher.matches("*46T>A") shouldBe true
            TransvarCDnaMatcher.matches("c.76_77delinsTT") shouldBe true
            TransvarCDnaMatcher.matches("76_77delinsTT") shouldBe true
        }

        "matches deletions" {
            TransvarCDnaMatcher.matches("c.76_78del") shouldBe true
            TransvarCDnaMatcher.matches("76_78del") shouldBe true
            TransvarCDnaMatcher.matches("c.76_78delACT") shouldBe true
            TransvarCDnaMatcher.matches("76_78delACT") shouldBe true
            TransvarCDnaMatcher.matches("c.301-3delT") shouldBe true
            TransvarCDnaMatcher.matches("301-3delT") shouldBe true
            TransvarCDnaMatcher.matches("c.5_9delAAGAG") shouldBe true
            TransvarCDnaMatcher.matches("5_9delAAGAG") shouldBe true
        }

        "matches duplications" {
            TransvarCDnaMatcher.matches("c.77_79dup") shouldBe true
            TransvarCDnaMatcher.matches("77_79dup") shouldBe true
            TransvarCDnaMatcher.matches("c.77_79dupCTG") shouldBe true
            TransvarCDnaMatcher.matches("77_79dupCTG") shouldBe true
        }

        "matches insertions" {
            TransvarCDnaMatcher.matches("c.76_77insT") shouldBe true
            TransvarCDnaMatcher.matches("76_77insT") shouldBe true
        }

        "matches indels" {
            TransvarCDnaMatcher.matches("c.112_117delinsTG") shouldBe true
            TransvarCDnaMatcher.matches("112_117delinsTG") shouldBe true
            TransvarCDnaMatcher.matches("c.112_117delAGGTCAinsTG") shouldBe true
            TransvarCDnaMatcher.matches("112_117delAGGTCAinsTG") shouldBe true
            TransvarCDnaMatcher.matches("c.113delinsTACTAGC") shouldBe true
            TransvarCDnaMatcher.matches("113delinsTACTAGC") shouldBe true
            TransvarCDnaMatcher.matches("c.113delGinsTACTAGC") shouldBe true
            TransvarCDnaMatcher.matches("113delGinsTACTAGC") shouldBe true
        }

        "does not match complex variants"{
            TransvarCDnaMatcher.matches("Fusion + c.324InsCGC") shouldBe false
            TransvarCDnaMatcher.matches("Fusion + 324InsCGC") shouldBe false
            TransvarCDnaMatcher.matches("Fusion + c.324InsCGC") shouldBe false
            TransvarCDnaMatcher.matches("Fusion + 324InsCGC") shouldBe false
        }

        "finds cdna variant in complex variant"{
            TransvarCDnaMatcher.contains("Splicing alteration (c.463+2C>T)") shouldBe true
            TransvarCDnaMatcher.contains("Splicing alteration (463+2C>T)") shouldBe true
            TransvarCDnaMatcher.contains("Fusion + c.324InsCGC") shouldBe true
            TransvarCDnaMatcher.contains("Fusion + 324InsCGC") shouldBe true
        }

        "does not match other variants" {
            TransvarCDnaMatcher.contains("DUX4 Fusions") shouldBe false
            TransvarCDnaMatcher.contains("SH2 DOMAIN MUTATION") shouldBe false
            TransvarCDnaMatcher.contains("DEL 485-490") shouldBe false
            TransvarCDnaMatcher.contains("RARE EX 18-21 MUT") shouldBe false
            TransvarCDnaMatcher.contains("Exon 12 splice site insertion") shouldBe false
            TransvarCDnaMatcher.contains("EWSR1-FLI1 Type 1") shouldBe false
            TransvarCDnaMatcher.contains("DPYD*2A HOMOZYGOSITY") shouldBe false
            TransvarCDnaMatcher.contains("3` UTR Polymorphism") shouldBe false
        }
    }
}
