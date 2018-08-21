package com.hartwig.hmftools.knowledgebaseimporter.transvar.matchers

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class TransvarGDnaMatcherTest : StringSpec() {

    init {
        "matches substitutions" {
            TransvarGDnaMatcher.matches("g.33038255C>A") shouldBe true
            TransvarGDnaMatcher.matches("g.33357783G>A") shouldBe true
        }

        "matches deletions" {
            TransvarGDnaMatcher.matches("g.19del") shouldBe true
            TransvarGDnaMatcher.matches("g.19_21del") shouldBe true
            TransvarGDnaMatcher.matches("g.19_21delTGC") shouldBe true
            TransvarGDnaMatcher.matches("g.32459297del") shouldBe true
        }

        "matches insertions" {
            TransvarGDnaMatcher.matches("g.32867861_32867862insT") shouldBe true
            TransvarGDnaMatcher.matches("g.32862923_32862924insCCT") shouldBe true
        }

        "matches duplications" {
            TransvarGDnaMatcher.matches("g.33229410dup") shouldBe true
            TransvarGDnaMatcher.matches("g.33229407_33229410dup") shouldBe true
            TransvarGDnaMatcher.matches("g.32862852_32862904dup") shouldBe true
            TransvarGDnaMatcher.matches("g.5dupT") shouldBe true
        }

        "matches indels" {
            TransvarGDnaMatcher.matches("g.6775delinsGA") shouldBe true
            TransvarGDnaMatcher.matches("g.6775_6777delinsC") shouldBe true
            TransvarGDnaMatcher.matches("g.9002_9009delinsTTT") shouldBe true
        }
    }
}
