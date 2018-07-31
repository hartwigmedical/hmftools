package com.hartwig.hmftools.idgenerator

import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class HmfIdTest : StringSpec({

    "patient index should be padded to 6 places" {
        HmfId("Any", 1).patientId shouldBe "HMF000001"
        HmfId("None", 999999).patientId shouldBe "HMF999999"
    }

    "large index is ok" {
        HmfId("Any", 1000000).patientId shouldBe "HMF1000000"
    }

})
