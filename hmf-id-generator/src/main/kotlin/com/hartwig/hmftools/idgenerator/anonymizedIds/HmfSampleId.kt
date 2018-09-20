package com.hartwig.hmftools.idgenerator.anonymizedIds

data class HmfSampleId(val hashId: HashId, val hmfPatientId: HmfPatientId) : AnonymizedId by hashId {
    val plaintext = hmfPatientId.plaintext + (64 + hashId.id).toChar()
}
