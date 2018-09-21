package com.hartwig.hmftools.idgenerator.anonymizedIds

import com.hartwig.hmftools.idgenerator.Hash

data class HmfPatientId(val hashId: HashId) : AnonymizedId by hashId {
    companion object {
        operator fun invoke(hashString: String, id: Int) = HmfPatientId(HashId(Hash(hashString), id))
    }

    fun updateHash(hash: Hash) = copy(hashId = hashId.copy(hash = hash))

    val plaintext = "HMF" + id.toString().padStart(6, '0')
}
