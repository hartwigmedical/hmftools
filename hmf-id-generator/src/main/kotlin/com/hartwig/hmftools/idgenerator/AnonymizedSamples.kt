package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.idgenerator.anonymizedIds.HmfSampleId
import com.hartwig.hmftools.idgenerator.ids.CanonicalPatientId
import com.hartwig.hmftools.idgenerator.ids.SampleId

class AnonymizedSamples private constructor(val password: String, val sampleIds: Collection<HmfSampleId>,
                                            private val samplesInput: SamplesInput) :
        Collection<HmfSampleId> by sampleIds {

    companion object {
        operator fun invoke(password: String, hmfSampleIds: Collection<HmfSampleId>, samplesInput: SamplesInput): AnonymizedSamples {
            val sampleIds = hmfSampleIds.toSet().sorted()
            return AnonymizedSamples(password, sampleIds, samplesInput)
        }
    }

    private val generator = IdGenerator(password)
    private val hmfSampleIdPerHash = sampleIds.groupBy { it.hash }
    private val hmfSampleIdPerPatientHash = sampleIds.groupBy { it.hmfPatientId.hash }

    operator fun get(sampleId: SampleId): HmfSampleId? {
        val hmfSampleIdsForHash = hmfSampleIdPerHash[generator.hash(sampleId.id)]
        hmfSampleIdsForHash ?: return null
        return if (hmfSampleIdsForHash.size == 1) hmfSampleIdsForHash.first()
        else hmfSampleIdsForHash.find { it.hmfPatientId.hash == generator.hash(samplesInput.canonicalId(sampleId.patientId).id) }
    }

    operator fun get(patientId: CanonicalPatientId): List<HmfSampleId> =
            hmfSampleIdPerPatientHash[generator.hash(patientId.id)] ?: emptyList()

    fun anonymizedSamplesMap(): Map<HmfSampleId, HmfSampleId> {
        val anonymizedPatients = AnonymizedPatients(password, sampleIds.map { it.hmfPatientId }.distinct(), samplesInput)
        return anonymizedPatients.anonymizedPatientMap().flatMap { (patientId, canonicalId) ->
            val patientIdSamples = hmfSampleIdPerPatientHash[patientId.hash] ?: emptyList()
            val canonicalIdSamples = hmfSampleIdPerPatientHash[canonicalId.hash]?.associateBy { it.hash } ?: emptyMap()
            patientIdSamples.map { Pair(it, canonicalIdSamples[it.hash]!!) }
        }.toMap()

    }
}
