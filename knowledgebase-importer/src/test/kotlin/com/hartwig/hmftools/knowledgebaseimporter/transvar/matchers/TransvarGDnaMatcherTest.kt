package com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers

import com.hartwig.hmftools.knowledgebaseimporter.knowledgebases.events.HgvsVariantType
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class TransvarGDnaMatcherTest : StringSpec() {

    init {
        "matches substitutions" {
            TransvarGDnaMatcher.matches("g.33038255C>A") shouldBe true
            TransvarGDnaMatcher.matches("33038255C>A") shouldBe true
            TransvarGDnaMatcher.type("g.33038255C>A") shouldBe HgvsVariantType.SUBSTITUTION
            TransvarGDnaMatcher.type("33038255C>A") shouldBe HgvsVariantType.SUBSTITUTION
        }

        "matches deletions" {
            TransvarGDnaMatcher.matches("g.19del") shouldBe true
            TransvarGDnaMatcher.matches("19del") shouldBe true
            TransvarGDnaMatcher.type("g.19del") shouldBe HgvsVariantType.DELETION
            TransvarGDnaMatcher.type("19del") shouldBe HgvsVariantType.DELETION
            TransvarGDnaMatcher.matches("g.19_21del") shouldBe true
            TransvarGDnaMatcher.matches("g.19_21delTGC") shouldBe true
            TransvarGDnaMatcher.matches("g.32459297del") shouldBe true
        }

        "matches range deletions" {
            TransvarGDnaMatcher.matches("g.19_21del") shouldBe true
            TransvarGDnaMatcher.matches("19_21del") shouldBe true
            TransvarGDnaMatcher.type("g.19_21del") shouldBe HgvsVariantType.DELETION
            TransvarGDnaMatcher.type("19_21del") shouldBe HgvsVariantType.DELETION
            TransvarGDnaMatcher.matches("g.19_21delTGC") shouldBe true
        }

        "matches range deletions with bases" {
            TransvarGDnaMatcher.matches("g.19_21delTGC") shouldBe true
            TransvarGDnaMatcher.matches("19_21delTGC") shouldBe true
            TransvarGDnaMatcher.type("g.19_21delTGC") shouldBe HgvsVariantType.DELETION
            TransvarGDnaMatcher.type("19_21delTGC") shouldBe HgvsVariantType.DELETION
        }

        "matches insertion" {
            TransvarGDnaMatcher.matches("g.32867861_32867862insT") shouldBe true
            TransvarGDnaMatcher.matches("32867861_32867862insT") shouldBe true
            TransvarGDnaMatcher.type("g.32867861_32867862insT") shouldBe HgvsVariantType.INSERTION
            TransvarGDnaMatcher.type("32867861_32867862insT") shouldBe HgvsVariantType.INSERTION
        }

        "matches insertion with multiple bases"{
            TransvarGDnaMatcher.matches("g.32862923_32862924insCCT") shouldBe true
            TransvarGDnaMatcher.matches("32862923_32862924insCCT") shouldBe true
            TransvarGDnaMatcher.type("g.32862923_32862924insCCT") shouldBe HgvsVariantType.INSERTION
            TransvarGDnaMatcher.type("32862923_32862924insCCT") shouldBe HgvsVariantType.INSERTION
        }

        "matches duplication" {
            TransvarGDnaMatcher.matches("g.33229410dup") shouldBe true
            TransvarGDnaMatcher.matches("33229410dup") shouldBe true
            TransvarGDnaMatcher.type("g.33229410dup") shouldBe HgvsVariantType.DUPLICATION
            TransvarGDnaMatcher.type("33229410dup") shouldBe HgvsVariantType.DUPLICATION
        }

        "matches duplication with base" {
            TransvarGDnaMatcher.matches("g.5dupT") shouldBe true
            TransvarGDnaMatcher.matches("5dupT") shouldBe true
            TransvarGDnaMatcher.type("g.5dupT") shouldBe HgvsVariantType.DUPLICATION
            TransvarGDnaMatcher.type("5dupT") shouldBe HgvsVariantType.DUPLICATION
        }

        "matches duplication with range" {
            TransvarGDnaMatcher.matches("g.33229407_33229410dup") shouldBe true
            TransvarGDnaMatcher.matches("33229407_33229410dup") shouldBe true
            TransvarGDnaMatcher.type("g.33229407_33229410dup") shouldBe HgvsVariantType.DUPLICATION
            TransvarGDnaMatcher.type("33229407_33229410dup") shouldBe HgvsVariantType.DUPLICATION
        }

        "matches indel" {
            TransvarGDnaMatcher.matches("g.6775delinsGA") shouldBe true
            TransvarGDnaMatcher.matches("6775delinsGA") shouldBe true
            TransvarGDnaMatcher.type("g.6775delinsGA") shouldBe HgvsVariantType.DELINS
            TransvarGDnaMatcher.type("6775delinsGA") shouldBe HgvsVariantType.DELINS
        }

        "matches indel with range" {
            TransvarGDnaMatcher.matches("g.6775_6777delinsC") shouldBe true
            TransvarGDnaMatcher.matches("6775_6777delinsC") shouldBe true
            TransvarGDnaMatcher.type("g.6775_6777delinsC") shouldBe HgvsVariantType.DELINS
            TransvarGDnaMatcher.type("6775_6777delinsC") shouldBe HgvsVariantType.DELINS
        }

        "matches indel with range and multiple bases" {
            TransvarGDnaMatcher.matches("g.9002_9009delinsTTT") shouldBe true
            TransvarGDnaMatcher.matches("9002_9009delinsTTT") shouldBe true
            TransvarGDnaMatcher.type("g.9002_9009delinsTTT") shouldBe HgvsVariantType.DELINS
            TransvarGDnaMatcher.type("9002_9009delinsTTT") shouldBe HgvsVariantType.DELINS
        }
    }
}
