package com.hartwig.hmftools.idgenerator

data class HmfId(val hash: String, val id: Int) {
    val patientId: String = "HMF" + id.toString().padStart(6, '0')
}