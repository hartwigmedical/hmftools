package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.idgenerator.anonymizedIds.HashId
import com.hartwig.hmftools.idgenerator.anonymizedIds.HmfPatientId
import com.hartwig.hmftools.idgenerator.anonymizedIds.HmfSampleId
import com.hartwig.hmftools.idgenerator.ids.PatientId
import com.hartwig.hmftools.idgenerator.ids.SampleId
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec

class AnonymizedSamplesTest : StringSpec() {
    private val password = "password"
    private val patient1 = PatientId("CPCT01990001")!!
    private val patient2 = PatientId("CPCT01990002")!!
    private val patient1Hash = Hash("6186c5ef806d9630ef209b8c32e6b7106cbc05e04bf39885b7d4663313217362")
    private val patient2Hash = Hash("c866a877e8bd5fe1b6a078953211369a190fdcd28a98934578beeb603f0d0125")
    private val sample1A = SampleId("${patient1.id}T")!!
    private val sample1B = SampleId("${patient1.id}TII")!!
    private val sample1AHash = Hash("2b838ee06b69f3880bca9d9675d0760aba8dede61e82492d1677e37fbcf20811")
    private val sample1BHash = Hash("453295f4e826bd4376f434a6f3d1fc151b30c5d8282af42f758f549fc962bb05")

    private val anonSample1A = HmfSampleId(HashId(sample1AHash, 1), HmfPatientId(HashId(patient1Hash, 1)))
    private val anonSample1B = HmfSampleId(HashId(sample1BHash, 2), HmfPatientId(HashId(patient1Hash, 1)))

    private val mappedAnonSample1A = HmfSampleId(HashId(sample1AHash, 1), HmfPatientId(HashId(patient2Hash, 1)))
    private val mappedAnonSample1B = HmfSampleId(HashId(sample1BHash, 2), HmfPatientId(HashId(patient2Hash, 1)))

    init {
        "resolves simple anonymized sample" {
            val samplesInput = SamplesInput(listOf(sample1A))
            val anonymizedSamples = AnonymizedSamples(password, listOf(anonSample1A), samplesInput)
            anonymizedSamples[sample1A] shouldBe anonSample1A
            anonymizedSamples[samplesInput.canonicalId(patient1)] shouldBe setOf(anonSample1A)
            anonymizedSamples.anonymizedSamplesMap().size shouldBe 0
        }

        "finds all anonymized samples for patient" {
            val samplesInput = SamplesInput(listOf(sample1A, sample1B))
            val anonymizedSamples = AnonymizedSamples(password, listOf(anonSample1A, anonSample1B), samplesInput)
            anonymizedSamples[samplesInput.canonicalId(patient1)] shouldBe setOf(anonSample1A, anonSample1B)
            anonymizedSamples.anonymizedSamplesMap().size shouldBe 0
        }

        "resolves single mapped anonymized sample" {
            val samplesInput = SamplesInput(listOf(sample1A), mapOf(patient1 to patient2))
            val anonymizedSamples = AnonymizedSamples(password, listOf(mappedAnonSample1A), samplesInput)
            anonymizedSamples[sample1A] shouldBe mappedAnonSample1A
            anonymizedSamples[samplesInput.canonicalId(patient1)] shouldBe setOf(mappedAnonSample1A)
            anonymizedSamples.anonymizedSamplesMap().size shouldBe 0
        }

        "resolves all anonymized samples for mapped patient" {
            val samplesInput = SamplesInput(listOf(sample1A, sample1B), mapOf(patient1 to patient2))
            val anonymizedSamples = AnonymizedSamples(password, listOf(mappedAnonSample1A, mappedAnonSample1B), samplesInput)
            anonymizedSamples[sample1A] shouldBe mappedAnonSample1A
            anonymizedSamples[sample1B] shouldBe mappedAnonSample1B
            anonymizedSamples[samplesInput.canonicalId(patient1)] shouldBe setOf(mappedAnonSample1A, mappedAnonSample1B)
            anonymizedSamples.anonymizedSamplesMap().size shouldBe 0
        }

        "resolves previously generatedId for removed sample" {
            val samplesInput = SamplesInput(listOf())
            val anonymizedSamples = AnonymizedSamples(password, listOf(anonSample1A), samplesInput)
            anonymizedSamples[sample1A] shouldBe anonSample1A
            anonymizedSamples[samplesInput.canonicalId(patient1)] shouldBe setOf(anonSample1A)
            anonymizedSamples.anonymizedSamplesMap().size shouldBe 0
        }

        "resolves correct sampleId for sample anonymized and mapped later" {
            val mappedAnonSample1A = mappedAnonSample1A.copy(hmfPatientId = HmfPatientId(HashId(patient2Hash, 2)))
            val samplesInput = SamplesInput(listOf(sample1A), mapOf(patient1 to patient2))
            val anonymizedSamples = AnonymizedSamples(password, listOf(anonSample1A, mappedAnonSample1A), samplesInput)
            anonymizedSamples[sample1A] shouldBe mappedAnonSample1A
            anonymizedSamples[samplesInput.canonicalId(patient1)] shouldBe setOf(mappedAnonSample1A)
            anonymizedSamples.anonymizedSamplesMap().size shouldBe 1
        }

        "resolves correct sampleId for patient anonymized, mapped and then deleted" {
            val mappedAnonSample1A = mappedAnonSample1A.copy(hmfPatientId = HmfPatientId(HashId(patient2Hash, 2)))
            val mappedAnonSample1B = mappedAnonSample1B.copy(hmfPatientId = HmfPatientId(HashId(patient2Hash, 2)))
            val samplesInput = SamplesInput(listOf(sample1B), mapOf(patient1 to patient2))
            val anonymizedSamples = AnonymizedSamples(password, listOf(anonSample1A, mappedAnonSample1A, mappedAnonSample1B), samplesInput)
            anonymizedSamples[sample1A] shouldBe mappedAnonSample1A
            anonymizedSamples[sample1B] shouldBe mappedAnonSample1B
            anonymizedSamples[samplesInput.canonicalId(patient1)] shouldBe setOf(mappedAnonSample1A, mappedAnonSample1B)
            anonymizedSamples.anonymizedSamplesMap().size shouldBe 1
        }

        "resolves original sampleId for patient anonymized, mapped and then unmapped" {
            val mappedAnonSample1A = mappedAnonSample1A.copy(hmfPatientId = HmfPatientId(HashId(patient2Hash, 2)))
            val samplesInput = SamplesInput(listOf(sample1A))
            val anonymizedSamples = AnonymizedSamples(password, listOf(anonSample1A, mappedAnonSample1A), samplesInput)
            anonymizedSamples[sample1A] shouldBe anonSample1A
            anonymizedSamples[samplesInput.canonicalId(patient1)] shouldBe setOf(anonSample1A)
            anonymizedSamples.anonymizedSamplesMap().size shouldBe 0
        }

        "resolves original sampleId for patient anonymized, mapped, unmapped and then deleted" {
            val mappedAnonSample1A = mappedAnonSample1A.copy(hmfPatientId = HmfPatientId(HashId(patient2Hash, 2)))
            val samplesInput = SamplesInput(listOf())
            val anonymizedSamples = AnonymizedSamples(password, listOf(anonSample1A, mappedAnonSample1A), samplesInput)
            anonymizedSamples[sample1A] shouldBe anonSample1A
            anonymizedSamples[samplesInput.canonicalId(patient1)] shouldBe setOf(anonSample1A)
            anonymizedSamples.anonymizedSamplesMap().size shouldBe 0
        }
    }
}
