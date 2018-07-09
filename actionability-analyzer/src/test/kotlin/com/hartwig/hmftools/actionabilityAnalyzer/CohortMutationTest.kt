package com.hartwig.hmftools.actionabilityAnalyzer

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class CohortMutationTest : StringSpec() {
    init {
        "splice oncogene is not potentially actionable" {
            val variant = CohortMutation("", "", "", "", "", "", "", "Splice", "",
                                         "SPLICE", "", "", "TRUE")
            variant.potentiallyActionable shouldBe false
        }

        "splice non-oncogene is potentially actionable" {
            val variant = CohortMutation("", "", "", "", "", "", "", "Splice", "",
                                         "SPLICE", "", "", "")
            variant.potentiallyActionable shouldBe true
        }

        "nonsense oncogene not potentially actionable" {
            val variant = CohortMutation("", "", "", "", "", "", "", "Nonsense", "",
                                         "NONSENSE_OR_FRAMESHIFT", "", "", "TRUE")
            variant.potentiallyActionable shouldBe false
        }

        "frameshift oncogene not potentially actionable" {
            val variant = CohortMutation("", "", "", "", "", "", "", "Frameshift", "",
                                         "NONSENSE_OR_FRAMESHIFT", "", "", "TRUE")
            variant.potentiallyActionable shouldBe false
        }

        "nonsense non-oncogene is potentially actionable" {
            val variant = CohortMutation("", "", "", "", "", "", "", "Nonsense", "",
                                         "NONSENSE_OR_FRAMESHIFT", "", "", "")
            variant.potentiallyActionable shouldBe true
        }

        "frameshift non-oncogene is potentially actionable" {
            val variant = CohortMutation("", "", "", "", "", "", "", "Frameshift", "",
                                         "NONSENSE_OR_FRAMESHIFT", "", "", "")
            variant.potentiallyActionable shouldBe true
        }

        "synonymous effect is not potentially actionable" {
            val variant = CohortMutation("", "", "", "", "", "", "", "Synonymous", "",
                                         "SYNONYMOUS", "", "", "")
            variant.potentiallyActionable shouldBe false
        }

        "missense is potentially actionable" {
            val variant = CohortMutation("", "", "", "", "", "", "", "Missense", "",
                                         "MISSENSE", "", "", "")
            variant.potentiallyActionable shouldBe true
        }

        "inframe is potentially actionable" {
            val variant = CohortMutation("", "", "", "", "", "", "", "Inframe", "",
                                         "", "", "", "")
            variant.potentiallyActionable shouldBe true
        }

        "missense oncogene is potentially actionable" {
            val variant = CohortMutation("", "", "", "", "", "", "", "Missense", "",
                                         "MISSENSE", "", "", "TRUE")
            variant.potentiallyActionable shouldBe true
        }
    }
}
