package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.idgenerator.ids.SampleId

data class AnonymizedRecord(val sampleId: String, val hmfSampleId: String) : CsvData {
    companion object {
        operator fun invoke(anonymizedSamples: AnonymizedSamples,  sampleId: SampleId): AnonymizedRecord {
            return AnonymizedRecord(sampleId.id, anonymizedSamples[sampleId]?.plaintext ?: "Unknown")
        }
    }
}
