package com.hartwig.hmftools.idgenerator.anonymizedIds

import com.hartwig.hmftools.idgenerator.Hash
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class HmfSampleIdTest : StringSpec() {
    init {
        "maps sample int id to letter"{
            HmfSampleId(1, 1, Hash("Any")).plaintext shouldBe "HMF000001A"
            HmfSampleId(1, 3, Hash("Any")).plaintext shouldBe "HMF000001C"
        }
    }
}
