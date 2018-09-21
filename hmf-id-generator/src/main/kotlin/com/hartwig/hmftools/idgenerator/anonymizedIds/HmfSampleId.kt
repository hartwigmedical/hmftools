package com.hartwig.hmftools.idgenerator.anonymizedIds

import com.hartwig.hmftools.idgenerator.Hash

data class HmfSampleId(val hashId: HashId, val hmfPatientId: HmfPatientId) : AnonymizedId by hashId {
    val plaintext = hmfPatientId.plaintext + (64 + hashId.id).toChar()
    fun updatePatientHash(hash: Hash) = copy(hmfPatientId = hmfPatientId.updateHash(hash))
    fun updateSampleHash(hash: Hash) = copy(hashId = hashId.copy(hash = hash))
}
