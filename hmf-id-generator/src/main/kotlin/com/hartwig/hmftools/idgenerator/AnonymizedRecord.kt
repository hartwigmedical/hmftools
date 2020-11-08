package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.common.amber.AmberAnonymous
import com.hartwig.hmftools.common.amber.ImmutableAmberAnonymous
import com.hartwig.hmftools.extensions.csv.CsvData
import com.hartwig.hmftools.idgenerator.anonymizedIds.HmfSampleId

data class AnonymizedRecord(val sampleId: String, val hmfSampleId: String) : CsvData {
    companion object {
        operator fun invoke(password: String, hmfSamples: List<HmfSampleId>, samples: List<String>): List<AnonymizedRecord> {
            val hmfSamplesMap: Map<Hash, HmfSampleId> = hmfSamples.map { Pair(it.hash, it) }.toMap()
            val generator = IdGenerator(password)
            return samples.map { x -> AnonymizedRecord(x, hmfSamplesMap[generator.hash(x)]!!.plaintext) }
        }
    }

    fun toAmberAnonymous(): AmberAnonymous = ImmutableAmberAnonymous.builder().sampleId(sampleId).hmfSampleId(hmfSampleId).deleted(false).build()
}
