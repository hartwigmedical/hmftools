package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.common.amber.AmberAnonymous
import com.hartwig.hmftools.common.amber.ImmutableAmberAnonymous
import com.hartwig.hmftools.extensions.csv.CsvData

data class HmfSample(val patientId: Int, val sampleId: Int, val deleted: Boolean, val sample: String) : Comparable<HmfSample> {

    val hmfSample = hmfSample()

    companion object {
        const val prefix = "HMF"

        operator fun invoke(amberAnonymous: AmberAnonymous): HmfSample {
            val hmfSampleId = amberAnonymous.hmfSampleId()

            val sampleId = hmfSampleId[hmfSampleId.length - 1].toInt() - 64
            val patientId = hmfSampleId.substring(prefix.length, hmfSampleId.length - 1).toInt()

            return HmfSample(patientId, sampleId, amberAnonymous.deleted(), amberAnonymous.sampleId())
        }
    }

    private fun hmfSample(): String {
        val patientString = prefix + patientId.toString().padStart(6, '0')
        val sampleString = (64 + sampleId).toChar()

        return patientString + sampleString.toString()
    }

    fun toCsv(generator: IdGenerator): HmfSampleCsv {
        return HmfSampleCsv(patientId.toString(), sampleId.toString(), hmfSample(), deleted.toString(), generator.hash(sample))
    }

    fun toAmberAnonymous(): AmberAnonymous {
        return ImmutableAmberAnonymous.builder().hmfSampleId(hmfSample()).sampleId(sample).deleted(deleted).build()
    }

    override fun compareTo(other: HmfSample): Int {
        return compareBy<HmfSample> { it.patientId }.thenComparing { hmfSampleId -> hmfSampleId.sampleId }.compare(this, other)
    }
}

data class HmfSampleCsv(val patientId: String, val sampleId: String, val hmfSampleId: String, val deleted: String, val sampleHash: String) : CsvData