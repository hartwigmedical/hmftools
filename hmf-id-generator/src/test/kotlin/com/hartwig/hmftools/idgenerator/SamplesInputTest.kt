package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.idgenerator.ids.PatientId
import com.hartwig.hmftools.idgenerator.ids.SampleId
import io.kotlintest.matchers.shouldBe
import io.kotlintest.matchers.shouldThrow
import io.kotlintest.specs.StringSpec

class SamplesInputTest : StringSpec() {
    private val patient1 = PatientId("CPCT01990001")!!
    private val patient2 = PatientId("CPCT01990002")!!
    private val patient3 = PatientId("CPCT01990003")!!
    private val sample1A = SampleId("${patient1.id}T")!!
    private val sample1B = SampleId("${patient1.id}TII")!!
    private val sample2A = SampleId("${patient2.id}T")!!
    private val sample3A = SampleId("${patient3.id}T")!!

    init {
        "does not allow chained renames" {
            val renames = mapOf(patient1 to patient2, patient2 to patient3)
            shouldThrow<IllegalStateException> { SamplesInput(listOf(sample1A, sample2A, sample3A), renames) }
        }

        "returns correct ids for a patient with no renames"{
            val renames = emptyMap<PatientId, PatientId>()
            val input = SamplesInput(listOf(sample1A, sample2A, sample3A), renames)
            input.patientIds(patient1) shouldBe setOf(patient1)
        }

        "returns correct ids for a patient with 2 ids"{
            val renames = mapOf(patient1 to patient2)
            val input = SamplesInput(listOf(sample1A, sample2A, sample3A), renames)
            input.patientIds(patient1) shouldBe setOf(patient1, patient2)
            input.patientIds(patient2) shouldBe input.patientIds(patient1)
        }

        "returns correct ids for a patient with 3 ids"{
            val renames = mapOf(patient1 to patient2, patient3 to patient2)
            val input = SamplesInput(listOf(sample1A, sample2A, sample3A), renames)
            input.patientIds(patient1) shouldBe setOf(patient1, patient2, patient3)
            input.patientIds(patient2) shouldBe input.patientIds(patient1)
            input.patientIds(patient3) shouldBe input.patientIds(patient1)
        }

        "returns correct samples for a patient with no renames" {
            val renames = emptyMap<PatientId, PatientId>()
            val input = SamplesInput(listOf(sample1A, sample1B, sample2A, sample3A), renames)
            input.sampleIds(patient1) shouldBe setOf(sample1A, sample1B)
        }

        "returns correct samples for a patient with 2 ids" {
            val renames = mapOf(patient1 to patient2)
            val input = SamplesInput(listOf(sample1A, sample1B, sample2A, sample3A), renames)
            input.sampleIds(patient1) shouldBe setOf(sample1A, sample1B, sample2A)
            input.sampleIds(patient2) shouldBe input.sampleIds(patient1)
        }

        "returns correct samples for a patient with 3 ids" {
            val renames = mapOf(patient1 to patient2, patient3 to patient2)
            val input = SamplesInput(listOf(sample1A, sample1B, sample2A, sample3A), renames)
            input.sampleIds(patient1) shouldBe setOf(sample1A, sample1B, sample2A, sample3A)
            input.sampleIds(patient2) shouldBe input.sampleIds(patient1)
            input.sampleIds(patient3) shouldBe input.sampleIds(patient1)
        }
    }
}
