package com.hartwig.hmftools.idgenerator.ids

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class SampleIdTest : StringSpec() {
    private val cpctPatient = "CPCT01990001"
    private val drupPatient = "DRUP01990002"
    private val cpctSample = "${cpctPatient}T"
    private val drupSample = "${drupPatient}T"

    init {
        "can create valid ids" {
            SampleId(cpctSample) shouldBe SampleId(PatientId(IdType.CPCT, cpctPatient), cpctSample)
            SampleId(drupSample) shouldBe SampleId(PatientId(IdType.DRUP, drupPatient), drupSample)
        }

        "ids that do not start with CPCT/DRUP are invalid" {
            SampleId("TEST01990001T") shouldBe null
        }

        "ids with incorrect patientId are invalid"{
            SampleId("CPCT0199001T") shouldBe null
            SampleId("CPCT0199001TIII") shouldBe null
            SampleId("CPCT019900010T") shouldBe null
        }
    }
}
