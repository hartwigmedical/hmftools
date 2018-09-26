package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.idgenerator.anonymizedIds.HashId
import com.hartwig.hmftools.idgenerator.anonymizedIds.HmfPatientId
import com.hartwig.hmftools.idgenerator.anonymizedIds.HmfSampleId

data class HmfSampleIdRecord(val patientId: String, val patientHash: String, val sampleId: String, val sampleHash: String,
                             val hmfSampleId: String) : CsvData {
    companion object {
        operator fun invoke(hmfSampleId: HmfSampleId): HmfSampleIdRecord {
            return HmfSampleIdRecord(hmfSampleId.hmfPatientId.id.toString(), hmfSampleId.hmfPatientId.hash.value, hmfSampleId.id.toString(),
                                     hmfSampleId.hash.value, hmfSampleId.plaintext)
        }
    }

    fun toHmfSampleId() = HmfSampleId(HashId(sampleHash, sampleId.toInt()), HmfPatientId(HashId(patientHash, patientId.toInt())))
}
