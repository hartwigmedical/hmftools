package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.idgenerator.anonymizedIds.HmfSampleId
import com.hartwig.hmftools.idgenerator.ids.PatientId
import com.hartwig.hmftools.idgenerator.ids.SampleId
import io.kotlintest.matchers.shouldBe
import io.kotlintest.matchers.shouldNotBe
import io.kotlintest.specs.StringSpec

private const val PASSWORD1 = "password"
private const val PASSWORD2 = "password_2"

class SampleAnonymizerTest : StringSpec() {
    private val anonymizer = SampleAnonymizer(PASSWORD1)

    private val patient1 = PatientId("CPCT01990001")!!
    private val patient2 = PatientId("CPCT01990002")!!
    private val patient3 = PatientId("CPCT01990003")!!
    private val sample1A = SampleId("${patient1.id}T")!!
    private val sample1B = SampleId("${patient1.id}TII")!!
    private val sample2A = SampleId("${patient2.id}T")!!
    private val sample3A = SampleId("${patient3.id}T")!!

    init {
        "generates sample id" {
            val samplesInput = SamplesInput(listOf(sample1A))
            val output = anonymizer.anonymize(PASSWORD1, samplesInput, emptyList())
            output.size shouldBe 1
            output[sample1A]!!.id shouldBe 1
        }

        "updates sample and patient hash on password change" {
            val samplesInput = SamplesInput(listOf(sample1A))
            val output = anonymizer.anonymize(PASSWORD1, samplesInput, emptyList())
            val updatedOutput = anonymizer.anonymize(PASSWORD2, samplesInput, output)
            updatedOutput.size shouldBe output.size

            val sampleId = output[sample1A]!!
            val updatedSampleId = updatedOutput[sample1A]!!
            updatedSampleId.id shouldBe sampleId.id
            updatedSampleId.hmfPatientId.id shouldBe sampleId.hmfPatientId.id
            updatedSampleId.hmfPatientId.hash shouldNotBe sampleId.hmfPatientId.hash
            collisions(updatedOutput, output) shouldBe 0
        }

        "multiple samples for same patient have different anonymized ids but same anonymized patientId" {
            val oneSampleInput = SamplesInput(listOf(sample1A))
            val twoSamplesInput = SamplesInput(listOf(sample1A, sample1B))
            val output = anonymizer.anonymize(PASSWORD1, oneSampleInput, emptyList())
            val updatedOutput = anonymizer.anonymize(PASSWORD1, twoSamplesInput, output)
            updatedOutput.size shouldBe output.size + 1

            val first = output[sample1A]!!
            val second = updatedOutput[sample1B]!!
            second.id shouldBe first.id + 1
            second.hash shouldNotBe first.hash
            second.hmfPatientId shouldBe first.hmfPatientId
        }

        "removing one already anonymized sample does not change the output" {
            val twoSamplesInput = SamplesInput(listOf(sample1A, sample1B))
            val oneSampleInput = SamplesInput(listOf(sample1A))
            val output = anonymizer.anonymize(PASSWORD1, twoSamplesInput, emptyList())
            val updatedOutput = anonymizer.anonymize(PASSWORD1, oneSampleInput, output)
            updatedOutput.sampleIds shouldBe output.sampleIds
            updatedOutput[sample1B]!!.id shouldBe 2
        }

        "changing password and removing sample for same patient updates patient hash but not sample hash for removed sample" {
            val twoSamplesInput = SamplesInput(listOf(sample1A, sample1B))
            val oneSampleInput = SamplesInput(listOf(sample1A))
            val output = anonymizer.anonymize(PASSWORD1, twoSamplesInput, emptyList())
            val updatedOutput = anonymizer.anonymize(PASSWORD2, oneSampleInput, output)
            updatedOutput.sampleIds.size shouldBe output.sampleIds.size
            collisions(updatedOutput, output) shouldBe 1
            updatedOutput.map { it.id } shouldBe output.map { it.id }
            updatedOutput.map { it.hmfPatientId.id }.toSet() shouldBe output.map { it.hmfPatientId.id }.toSet()
            updatedOutput.map { it.hmfPatientId.hash }.toSet() shouldNotBe output.map { it.hmfPatientId.hash }.toSet()
        }

        "replacing an already anonymized sample adds entry to the output" {
            val twoSamplesInput = SamplesInput(listOf(sample1A, sample1B))
            val replacedSampleInput = SamplesInput(listOf(sample1A, sample2A))
            val output = anonymizer.anonymize(PASSWORD1, twoSamplesInput, emptyList())
            val updatedOutput = anonymizer.anonymize(PASSWORD1, replacedSampleInput, output)
            updatedOutput.sampleIds.size shouldBe output.sampleIds.size + 1
            collisions(updatedOutput, output) shouldBe 2
        }

        "samples for patients mapped from start generate single hmfPatientId" {
            val patientMap = mapOf(patient1 to patient2)
            val samplesInput = SamplesInput(listOf(sample1A, sample2A), patientMap)
            val output = anonymizer.anonymize(PASSWORD1, samplesInput, emptyList())
            output.size shouldBe 2
            output[sample1A]!!.id shouldBe 1
            output[sample2A]!!.id shouldBe 2
            output[sample1A]!!.hmfPatientId shouldBe output[sample2A]!!.hmfPatientId
            output.anonymizedSamplesMap().size shouldBe 0
        }

        "samples for patients mapped from start generate single hmfPatientId after password change" {
            val patientMap = mapOf(patient1 to patient2)
            val samplesInput = SamplesInput(listOf(sample1A, sample2A), patientMap)
            val output = anonymizer.anonymize(PASSWORD1, samplesInput, emptyList())
            val updatedOutput = anonymizer.anonymize(PASSWORD2, samplesInput, output)

            updatedOutput.size shouldBe 2
            updatedOutput[sample1A]!!.id shouldBe 1
            updatedOutput[sample2A]!!.id shouldBe 2
            updatedOutput[sample1A]!!.hmfPatientId shouldBe updatedOutput[sample2A]!!.hmfPatientId
            updatedOutput[sample1A]!!.hmfPatientId.id shouldBe output[sample1A]!!.hmfPatientId.id
            updatedOutput[sample1A]!!.hmfPatientId.hash shouldNotBe output[sample1A]!!.hmfPatientId.hash
            collisions(updatedOutput, output) shouldBe 0
            output.anonymizedSamplesMap().size shouldBe 0
        }

        "samples for patients mapped later add an extra item to output" {
            val patientMap = mapOf(patient1 to patient2)
            val initialInput = SamplesInput(listOf(sample1A, sample2A))
            val mappedInput = SamplesInput(listOf(sample1A, sample2A), patientMap)
            val output = anonymizer.anonymize(PASSWORD1, initialInput, emptyList())
            output.size shouldBe 2

            val updatedOutput = anonymizer.anonymize(PASSWORD1, mappedInput, output)
            updatedOutput.size shouldBe 3
            updatedOutput[sample1A]!!.id shouldBe 2
            updatedOutput[sample2A]!!.id shouldBe 1
            updatedOutput[sample1A]!!.hmfPatientId shouldBe updatedOutput[sample2A]!!.hmfPatientId
            updatedOutput[sample1A]!!.hash shouldBe output[sample1A]!!.hash
            updatedOutput[sample2A]!!.hash shouldBe output[sample2A]!!.hash
            updatedOutput[sample1A]!!.hmfPatientId shouldNotBe output[sample1A]!!.hmfPatientId
            updatedOutput.anonymizedSamplesMap().size shouldBe 1
            updatedOutput.anonymizedSamplesMap().map { it.key.plaintext to it.value.plaintext }.toMap() shouldBe mapOf("HMF000001A" to "HMF000002B")
        }

        "samples for patients mapped later correctly update hashes on password change" {
            val patientMap = mapOf(patient1 to patient2)
            val initialInput = SamplesInput(listOf(sample1A, sample2A))
            val mappedInput = SamplesInput(listOf(sample1A, sample2A), patientMap)
            val output = anonymizer.anonymize(PASSWORD1, initialInput, emptyList())
            output.size shouldBe 2

            val updatedOutput = anonymizer.anonymize(PASSWORD2, mappedInput, output)
            updatedOutput.size shouldBe 3
            updatedOutput[sample1A]!!.id shouldBe 2
            updatedOutput[sample2A]!!.id shouldBe 1
            updatedOutput[sample1A]!!.hmfPatientId shouldBe updatedOutput[sample2A]!!.hmfPatientId
            updatedOutput[sample1A]!!.hash shouldNotBe output[sample1A]!!.hash
            updatedOutput[sample2A]!!.hash shouldNotBe output[sample2A]!!.hash
            updatedOutput[sample1A]!!.hmfPatientId shouldNotBe output[sample1A]!!.hmfPatientId
            updatedOutput.map { it.hmfPatientId.hash }.none { output.map { it.hmfPatientId.hash }.contains(it) } shouldBe true
            collisions(updatedOutput, output) shouldBe 0
            updatedOutput.anonymizedSamplesMap().size shouldBe 1
            updatedOutput.anonymizedSamplesMap().map { it.key.plaintext to it.value.plaintext }.toMap() shouldBe mapOf("HMF000001A" to "HMF000002B")
        }

        "samples added for patients mapped later are only added with the canonical hmfPatientId" {
            val patientMap = mapOf(patient1 to patient2)
            val initialInput = SamplesInput(listOf(sample1A, sample2A))
            val mappedInput = SamplesInput(listOf(sample1A, sample2A, sample1B), patientMap)
            val output = anonymizer.anonymize(PASSWORD1, initialInput, emptyList())
            output.size shouldBe 2

            val updatedOutput = anonymizer.anonymize(PASSWORD1, mappedInput, output)
            updatedOutput.size shouldBe 4
            updatedOutput[sample1A]!!.id shouldBe 2
            updatedOutput[sample1B]!!.id shouldBe 3
            updatedOutput[sample2A]!!.id shouldBe 1
            updatedOutput[sample1A]!!.hmfPatientId shouldBe updatedOutput[sample2A]!!.hmfPatientId
            updatedOutput[sample1B]!!.hmfPatientId shouldBe updatedOutput[sample2A]!!.hmfPatientId
            updatedOutput[sample1A]!!.hmfPatientId shouldNotBe output[sample1A]!!.hmfPatientId
            updatedOutput.anonymizedSamplesMap().size shouldBe 1
            updatedOutput.anonymizedSamplesMap().map { it.key.plaintext to it.value.plaintext }.toMap() shouldBe mapOf("HMF000001A" to "HMF000002B")
        }

        "samples added for patients mapped later that are then unmapped continue where they left off" {
            val patientMap = mapOf(patient1 to patient2)
            val initialInput = SamplesInput(listOf(sample1A, sample2A))
            val mappedInput = SamplesInput(listOf(sample1A, sample2A, sample1B), patientMap)
            val unmappedInput = SamplesInput(listOf(sample1A, sample2A, sample1B))
            val output = anonymizer.anonymize(PASSWORD1, initialInput, emptyList())
            output.size shouldBe 2
            val output2 = anonymizer.anonymize(PASSWORD1, mappedInput, output)
            output2.size shouldBe 4
            val output3 = anonymizer.anonymize(PASSWORD1, unmappedInput, output2)
            output3.size shouldBe 5
            output3[sample1B]!!.id shouldBe 2
            output3[sample1B]!!.hmfPatientId shouldBe output[sample1A]!!.hmfPatientId
            output3.anonymizedSamplesMap().size shouldBe 0
        }

        "samples added for patients mapped later that are then unmapped continue where they left off with password changes" {
            val patientMap = mapOf(patient1 to patient2)
            val initialInput = SamplesInput(listOf(sample1A, sample2A))
            val mappedInput = SamplesInput(listOf(sample1A, sample2A, sample1B), patientMap)
            val unmappedInput = SamplesInput(listOf(sample1A, sample2A, sample1B))
            val output = anonymizer.anonymize(PASSWORD1, initialInput, emptyList())
            output.size shouldBe 2
            val output2 = anonymizer.anonymize(PASSWORD2, mappedInput, output)
            output2.size shouldBe 4
            val newAnonymizer = SampleAnonymizer(PASSWORD2)
            val output3 = newAnonymizer.anonymize("pass3", unmappedInput, output2)
            output3.size shouldBe 5
            output3[sample1B]!!.id shouldBe 2
            output3[sample1B]!!.hmfPatientId.id shouldBe output[sample1A]!!.hmfPatientId.id
            collisions(output2, output) shouldBe 0
            collisions(output3, output2) shouldBe 0
            collisions(output3, output) shouldBe 0
            output3.anonymizedSamplesMap().size shouldBe 0
        }

        "samples for patients mapped later that are then removed do not end up in anonymized samples map" {
            val patientMap = mapOf(patient1 to patient2)
            val initialInput = SamplesInput(listOf(sample1A, sample2A))
            val mappedInput = SamplesInput(listOf(sample1A, sample2A, sample1B), patientMap)
            val mappedInputWithRemovedSample = SamplesInput(listOf(sample1A, sample2A), patientMap)
            val output = anonymizer.anonymize(PASSWORD1, initialInput, emptyList())
            output.size shouldBe 2
            val output2 = anonymizer.anonymize(PASSWORD1, mappedInput, output)
            output2.size shouldBe 4
            val output3 = anonymizer.anonymize(PASSWORD1, mappedInputWithRemovedSample, output2)
            output3.size shouldBe 4
            output3[sample1B]!!.id shouldBe 3
            output3[sample1B]!!.hmfPatientId shouldBe output[sample2A]!!.hmfPatientId
            output3.anonymizedSamplesMap().size shouldBe 1
        }

        "samples for patients mapped later that are then unmapped and samples removed do not end up in anonymized samples map" {
            val patientMap = mapOf(patient1 to patient2)
            val initialInput = SamplesInput(listOf(sample1A, sample2A))
            val mappedInput = SamplesInput(listOf(sample1A, sample2A, sample1B), patientMap)
            val unmappedInputWithRemovedSample = SamplesInput(listOf(sample1A, sample2A))
            val output = anonymizer.anonymize(PASSWORD1, initialInput, emptyList())
            output.size shouldBe 2
            val output2 = anonymizer.anonymize(PASSWORD1, mappedInput, output)
            output2.size shouldBe 4
            val output3 = anonymizer.anonymize(PASSWORD1, unmappedInputWithRemovedSample, output2)
            output3.size shouldBe 4
            output3[sample1B]!!.plaintext shouldBe "HMF000002C"
            output3.anonymizedSamplesMap().size shouldBe 0
        }

        "samples for patient mapped later and then unmapped resolve to the original id" {
            val patientMap = mapOf(patient1 to patient2)
            val initialInput = SamplesInput(listOf(sample1A, sample2A))
            val mappedInput = SamplesInput(listOf(sample1A, sample2A), patientMap)
            val unmappedInput = SamplesInput(listOf(sample1A, sample2A))
            val output = anonymizer.anonymize(PASSWORD1, initialInput, emptyList())
            output.size shouldBe 2
            output[sample1A]!!.plaintext shouldBe "HMF000001A"
            val output2 = anonymizer.anonymize(PASSWORD1, mappedInput, output)
            output2.size shouldBe 3
            output2[sample1A]!!.plaintext shouldBe "HMF000002B"
            val output3 = anonymizer.anonymize(PASSWORD1, unmappedInput, output2)
            output3.size shouldBe 3
            output3[sample1A]!!.plaintext shouldBe "HMF000001A"
        }
    }

    private fun collisions(samples: Collection<HmfSampleId>, otherSamples: Collection<HmfSampleId>): Int {
        val patientHashes = samples.map { it.hash }.toSet()
        val otherHashes = otherSamples.map { it.hash }.toSet()
        return patientHashes.count { otherHashes.contains(it) }
    }
}
