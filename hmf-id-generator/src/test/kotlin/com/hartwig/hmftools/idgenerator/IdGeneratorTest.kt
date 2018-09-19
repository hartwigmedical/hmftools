package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.idgenerator.anonymizedIds.HashId
import io.kotlintest.matchers.shouldBe
import io.kotlintest.properties.forAll
import io.kotlintest.properties.forNone
import io.kotlintest.specs.StringSpec

private const val PASSWORD1 = "password"
private const val PASSWORD2 = "password_2"

class IdGeneratorTest : StringSpec() {
    private val generator = IdGenerator(PASSWORD1)
    private val generator2 = IdGenerator(PASSWORD2)

    init {
        "same password, same input generates same ids" {
            forAll { a: String ->
                generator.generate(listOf(a)) == generator.generate(listOf(a))
            }
        }

        "different passwords generate different ids" {
            forNone { a: String ->
                generator.generate(listOf(a)) == generator2.generate(listOf(a))
            }
        }

        "repeated samples get one id" {
            forAll { a: String, b: String ->
                val samples = listOf(a, b, b, a)
                val ids = generator.generate(samples)
                if (a != b) {
                    (ids.size == 2) and (hmfId(a, generator, ids) == 1) and (hmfId(b, generator, ids) == 2)
                } else {
                    (ids.size == 1) and (hmfId(a, generator, ids) == 1)
                }
            }
        }

        "different strings generate different ids, same strings generate same ids" {
            forAll { a: String, b: String ->
                if (a == b) {
                    generator.generate(listOf(a)) == generator.generate(listOf(b))
                } else {
                    generator.generate(listOf(a)) != generator.generate(listOf(b))
                }
            }
        }

        "ids are stable" {
            forAll { a: String, b: String, c: String ->
                val samples = listOf(a, b, c)
                val oldIds = generator.generate(samples)
                val newIds = generator.update(PASSWORD2, samples, oldIds)
                samples.all {
                    (hmfId(it, generator, oldIds) != null) and (hmfId(it, generator, oldIds) == hmfId(it, generator2, newIds))
                }
            }
        }

        "ids for new patients are incremental" {
            forAll { a: String, b: String, c: String ->
                val samples = listOf(a, b)
                val newSamples = samples + c
                val oldIds = generator.generate(samples)
                val newIds = generator.update(PASSWORD2, newSamples, oldIds)
                val previousMaxId = oldIds.maxBy { it.id }!!.id
                if (samples.contains(c)) {
                    if (c == a) {
                        hmfId(c, generator2, newIds)!! == 1
                    } else {
                        hmfId(c, generator2, newIds)!! == 2
                    }
                } else {
                    hmfId(c, generator2, newIds) == (previousMaxId + 1)
                }
            }
        }

        "can update when existing sample list is empty" {
            val newSample = "sample_A"
            val samples = emptyList<String>()
            val newSamples = samples + newSample
            val oldIds = generator.generate(samples)
            val newIds = generator.update(PASSWORD2, newSamples, oldIds)
            hmfId(newSample, generator2, newIds)!! shouldBe 1
        }

        "can update when existing sample list contains only one id" {
            val sampleA = "sample_A"
            val newSample = "sample_B"
            val samples = listOf(sampleA)
            val newSamples = samples + newSample
            val oldIds = generator.generate(samples)
            val newIds = generator.update(PASSWORD2, newSamples, oldIds)
            hmfId(newSample, generator2, newIds)!! shouldBe 2
        }

        "does not delete ids" {
            val samples = listOf("a", "b", "c")
            val newSamples = listOf("a", "c")
            val oldIds = generator.generate(samples).toSet()
            val newIds = generator.updateIdList(newSamples, oldIds).toSet()
            newIds shouldBe oldIds
        }

        "changing password does not delete ids" {
            val samples = listOf("a", "b", "c")
            val newSamples = listOf("a", "c")
            val oldIds = generator.generate(samples).toSet()
            val newIds = generator.changePassword(PASSWORD2, newSamples, oldIds).toSet()
            newIds.size shouldBe oldIds.size
        }

        "updating an empty id list is the same as generating new ids" {
            val samples = listOf("a", "b", "c")
            val oldIds = generator.generate(samples).toSet()
            val newIds = generator.updateIdList(samples, emptyList()).toSet()
            newIds shouldBe oldIds
        }

        "once an id gets assigned it will be reused" {
            val samples = listOf("a", "b", "c")
            val newSamples = listOf("a")
            val sampleIds = generator.generate(samples)
            val newSampleIds = generator.updateIdList(newSamples, sampleIds)
            val reUpdatedSampleIds = generator.updateIdList(samples, newSampleIds)
            reUpdatedSampleIds shouldBe sampleIds
        }

        "samples that are deleted before password change then re-added will get new id" {
            val deletedSample = "b"
            val samples = listOf("a", deletedSample, "c")
            val newSamples = samples - deletedSample
            val anonSamples = generator.generate(samples)
            val anonNewSamples = generator.changePassword(PASSWORD2, newSamples, anonSamples)
            val reAddedAnonSamples = generator2.updateIdList(samples, anonNewSamples)
            reAddedAnonSamples.size shouldBe (anonSamples.size + 1)
            hmfId(deletedSample, generator2, reAddedAnonSamples) shouldBe 4
        }

        "duplicated input generates single id" {
            val samples = listOf("a", "b")
            val newSamples = samples + "a"
            val anonSamples = generator.generate(samples)
            val anonNewSamples = generator.update(PASSWORD2, newSamples, anonSamples)
            anonNewSamples.size shouldBe newSamples.distinct().size
        }
    }

    private fun hmfId(sample: String, generator: IdGenerator, ids: Collection<HashId>): Int? {
        val sampleHash = generator.generate(listOf(sample)).first().hash
        return ids.associateBy { it.hash }[sampleHash]?.id
    }
}
