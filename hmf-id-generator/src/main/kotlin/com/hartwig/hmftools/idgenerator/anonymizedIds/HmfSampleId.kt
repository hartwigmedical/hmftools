package com.hartwig.hmftools.idgenerator.anonymizedIds

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.idgenerator.Hash

data class HmfSampleId(val patientId: Int, val sampleId: Int, val hash: Hash) : Comparable<HmfSampleId> {
    private val hmfPatientId = HmfPatientId(patientId)
    val plaintext = hmfPatientId.plaintext + (64 + sampleId).toChar()

    fun updateHash(hash: Hash) = copy(hash = hash)

    override fun compareTo(other: HmfSampleId): Int {
        return compareBy<HmfSampleId> { it.patientId }.thenComparing { hmfSampleId -> hmfSampleId.sampleId }.compare(this, other)
    }

    fun toCsv(): HmfSampleIdCsv {
        return HmfSampleIdCsv(patientId.toString(), sampleId.toString(), hash.value, plaintext)
    }
}

data class HmfSampleIdCsv(val patientId: String, val sampleId: String, val sampleHash: String, val hmfSampleId: String) : CsvData {
    fun toHmfSampleId(): HmfSampleId {
        return HmfSampleId(patientId.toInt(), sampleId.toInt(), Hash(sampleHash))
    }
}
