package com.hartwig.hmftools.idgenerator.anonymizedIds

import com.hartwig.hmftools.idgenerator.Hash

data class HmfSampleId(val hashId: HashId, val hmfPatientId: HmfPatientId) : AnonymizedId by hashId, Comparable<HmfSampleId> {
    val plaintext = hmfPatientId.plaintext + (64 + hashId.id).toChar()
    fun updatePatientHash(hash: Hash) = copy(hmfPatientId = hmfPatientId.updateHash(hash))
    fun updateSampleHash(hash: Hash) = copy(hashId = hashId.copy(hash = hash))

    override fun compareTo(other: HmfSampleId): Int {
        return compareBy<HmfSampleId> { it.hmfPatientId }.thenComparing { hmfSampleId -> hmfSampleId.hashId.id }.compare(this, other)
    }

    fun simplify(): HmfSampleIdSimple {
        return HmfSampleIdSimple(hmfPatientId.id, id, hash)
    }
}

data class HmfSampleIdSimple(val patientId: Int, val sampleId: Int, val hash: Hash): Comparable<HmfSampleIdSimple> {
    private val hmfPatientId = HmfPatientIdSimple(patientId)
    val plaintext = hmfPatientId.plaintext + (64 + sampleId).toChar()

    fun updateHash(hash: Hash) = copy(hash = hash)

    fun complicate(): HmfSampleId {
        return HmfSampleId(HashId(hash, sampleId), HmfPatientId("na", patientId))
    }

    override fun compareTo(other: HmfSampleIdSimple): Int {
        return compareBy<HmfSampleIdSimple> { it.patientId }.thenComparing { hmfSampleId -> hmfSampleId.sampleId }.compare(this, other)
    }
}
