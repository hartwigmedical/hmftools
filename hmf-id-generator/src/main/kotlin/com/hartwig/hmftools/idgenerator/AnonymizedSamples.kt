package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.idgenerator.anonymizedIds.HmfSampleId
import com.hartwig.hmftools.idgenerator.ids.CanonicalPatientId
import com.hartwig.hmftools.idgenerator.ids.SampleId

class AnonymizedSamples(val password: String, private val hmfSampleIds: Collection<HmfSampleId>, private val samplesInput: SamplesInput) :
        Collection<HmfSampleId> by hmfSampleIds {
    private val generator = IdGenerator(password)
    val sampleIds = hmfSampleIds.toSet().sorted()
    private val hmfSampleIdPerHash = hmfSampleIds.groupBy { it.hash }
    private val hmfSampleIdPerPatientHash = hmfSampleIds.groupBy { it.hmfPatientId.hash }

    operator fun get(sampleId: SampleId): HmfSampleId? = hmfSampleIdPerHash[generator.hash(sampleId.id)]
            ?.find { it.hmfPatientId.hash == generator.hash(samplesInput.canonicalId(sampleId.patientId).patientId.id) }

    operator fun get(patientId: CanonicalPatientId): List<HmfSampleId> =
            hmfSampleIdPerPatientHash[generator.hash(patientId.patientId.id)] ?: emptyList()
}
