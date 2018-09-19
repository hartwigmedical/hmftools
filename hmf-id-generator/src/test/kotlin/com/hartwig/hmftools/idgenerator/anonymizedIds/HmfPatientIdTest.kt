package com.hartwig.hmftools.idgenerator.anonymizedIds

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class HmfPatientIdTest : StringSpec() {

    init {
        "patient index should be padded to 6 places"{
            HmfPatientId("Any", 1).plaintext shouldBe "HMF000001"
            HmfPatientId("None", 999999).plaintext shouldBe "HMF999999"
        }

        "large index is ok"{
            HmfPatientId("Any", 1000000).plaintext shouldBe "HMF1000000"
        }
    }
}
