package com.hartwig.hmftools.knowledgebaseimporter.civic

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class CivicEvidenceTest : StringSpec() {

    init {
        "handles empty list" {
            CivicEvidence.splitDrugs("") shouldBe emptyList<String>()
        }

        "extracts single drug" {
            CivicEvidence.splitDrugs("Cetuximab") shouldBe listOf("Cetuximab")
        }

        "extracts single drug with comma" {
            CivicEvidence.splitDrugs("Cetuximab,") shouldBe listOf("Cetuximab")
        }

        "extracts single drug with whitespace and comma" {
            CivicEvidence.splitDrugs("Cetuximab  ,") shouldBe listOf("Cetuximab")
        }

        "splits simple drug list" {
            CivicEvidence.splitDrugs("Cetuximab, Panitumumab") shouldBe listOf("Cetuximab", "Panitumumab")
        }

        "does not split drug with examples"{
            CivicEvidence.splitDrugs("BEZ235 (NVP-BEZ235, Dactolisib)") shouldBe listOf("BEZ235 (NVP-BEZ235, Dactolisib)")
        }

        "correctly splits combination treatment" {
            CivicEvidence.splitDrugs("BEZ235 (NVP-BEZ235, Dactolisib), Cetuximab") shouldBe
                    listOf("BEZ235 (NVP-BEZ235, Dactolisib)", "Cetuximab")
        }
    }
}
