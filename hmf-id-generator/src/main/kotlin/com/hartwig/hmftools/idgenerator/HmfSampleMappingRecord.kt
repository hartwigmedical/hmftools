package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.idgenerator.anonymizedIds.HmfSampleId

data class HmfSampleMappingRecord(val mappedHmfSampleId: String, val hmfSampleId: String) : CsvData {
    companion object {
        operator fun invoke(mappedHmfSampleId: HmfSampleId, hmfSampleId: HmfSampleId): HmfSampleMappingRecord {
            return HmfSampleMappingRecord(mappedHmfSampleId.plaintext, hmfSampleId.plaintext)
        }
    }
}
