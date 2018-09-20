package com.hartwig.hmftools.idgenerator.ids

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class PatientIdTest : StringSpec() {
    private val cpctId = "CPCT01990001"
    private val drupId = "DRUP01990002"

    init {
        "can create valid ids" {
            PatientId(cpctId) shouldBe PatientId(IdType.CPCT, cpctId)
            PatientId(drupId) shouldBe PatientId(IdType.DRUP, drupId)
        }

        "ids that do not start with CPCT/DRUP are invalid" {
            PatientId("TEST01990001") shouldBe null
        }

        "ids with length <> 12 are invalid" {
            PatientId("CPCT0199001") shouldBe null
            PatientId("CPCT019900010") shouldBe null
        }
    }
}
