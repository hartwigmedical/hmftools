package com.hartwig.hmftools.idgenerator.anonymizedIds

data class HmfPatientId(val hashId: HashId) : AnonymizedId by hashId {
    val plaintext = "HMF" + id.toString().padStart(6, '0')
}
