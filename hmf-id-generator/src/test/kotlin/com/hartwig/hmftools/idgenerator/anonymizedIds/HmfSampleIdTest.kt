package com.hartwig.hmftools.idgenerator.anonymizedIds

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class HmfSampleIdTest : StringSpec() {
    private val patientId = HmfPatientId("Any", 1)

    init {
        "maps sample int id to letter"{
            HmfSampleId(HashId("Any", 1), patientId).plaintext shouldBe "HMF000001A"
            HmfSampleId(HashId("Any", 3), patientId).plaintext shouldBe "HMF000001C"
        }
    }
}
