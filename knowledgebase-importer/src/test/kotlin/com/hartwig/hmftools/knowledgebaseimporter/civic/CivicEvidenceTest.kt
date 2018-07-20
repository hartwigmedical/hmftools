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

        "extracts single drug with space" {
            CivicEvidence.splitDrugs("Cetuximab phosphate") shouldBe listOf("Cetuximab phosphate")
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

        "does not split complicated drug with both round and square brackets"{
            CivicEvidence.splitDrugs("Bis-2-[5-(phenylacetamide)-1,3,4-thiadiazol-2-yl]ethyl Sulfide") shouldBe
                    listOf("Bis-2-[5-(phenylacetamide)-1,3,4-thiadiazol-2-yl]ethyl Sulfide")

        }

        "correctly splits combination treatment" {
            CivicEvidence.splitDrugs("BEZ235 (NVP-BEZ235, Dactolisib), Cetuximab") shouldBe
                    listOf("BEZ235 (NVP-BEZ235, Dactolisib)", "Cetuximab")
        }

        "correctly splits combination treatment with complicated drug" {
            CivicEvidence.splitDrugs("Bis-2-[5-(phenylacetamide)-1,3,4-thiadiazol-2-yl]ethyl Sulfide, BEZ235 (NVP-BEZ235, Dactolisib), Cetuximab") shouldBe
                    listOf("Bis-2-[5-(phenylacetamide)-1,3,4-thiadiazol-2-yl]ethyl Sulfide", "BEZ235 (NVP-BEZ235, Dactolisib)", "Cetuximab")
        }
    }
}
