package com.hartwig.hmftools.actionabilityAnalyzer

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class CohortMutationTest : StringSpec() {
    init {
        "splice oncogene is not potentially actionable" {
            val variant = CohortMutation("", "", "", "", "", "", "", "", "",
                                         "SPLICE", "", "TRUE")
            variant.potentiallyActionable shouldBe false
        }

        "splice non-oncogene is potentially actionable" {
            val variant = CohortMutation("", "", "", "", "", "", "", "", "",
                                         "SPLICE", "", "")
            variant.potentiallyActionable shouldBe true
        }

        "nonsense or frameshift oncogene not potentially actionable" {
            val variant = CohortMutation("", "", "", "", "", "", "", "", "",
                                         "NONSENSE_OR_FRAMESHIFT", "", "TRUE")
            variant.potentiallyActionable shouldBe false
        }

        "nonsense or frameshift non-oncogene is potentially actionable" {
            val variant = CohortMutation("", "", "", "", "", "", "", "", "",
                                         "NONSENSE_OR_FRAMESHIFT", "", "")
            variant.potentiallyActionable shouldBe true
        }

        "synonymous effect is not potentially actionable" {
            val variant = CohortMutation("", "", "", "", "", "", "", "", "",
                                         "SYNONYMOUS", "", "")
            variant.potentiallyActionable shouldBe false
        }

        "none effect is not potentially actionable" {
            val variant = CohortMutation("", "", "", "", "", "", "", "", "",
                                         "NONE", "", "")
            variant.potentiallyActionable shouldBe false
        }

        "missense is potentially actionable" {
            val variant = CohortMutation("", "", "", "", "", "", "", "", "",
                                         "MISSENSE", "", "")
            variant.potentiallyActionable shouldBe true
        }

        "missense oncogene is potentially actionable" {
            val variant = CohortMutation("", "", "", "", "", "", "", "", "",
                                         "MISSENSE", "", "TRUE")
            variant.potentiallyActionable shouldBe true
        }
    }
}