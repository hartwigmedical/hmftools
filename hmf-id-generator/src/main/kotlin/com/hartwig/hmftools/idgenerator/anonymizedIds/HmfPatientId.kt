package com.hartwig.hmftools.idgenerator.anonymizedIds

data class HmfPatientId(val id: Int) : Comparable<HmfPatientId> {
    val plaintext = "HMF" + id.toString().padStart(6, '0')

    override fun compareTo(other: HmfPatientId): Int {
        return id.compareTo(other.id)
    }
}
