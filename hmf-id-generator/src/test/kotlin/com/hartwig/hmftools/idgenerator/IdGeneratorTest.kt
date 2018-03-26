package com.hartwig.hmftools.idgenerator

import io.kotlintest.properties.forAll
import io.kotlintest.properties.forNone
import io.kotlintest.specs.StringSpec

const val PASSWORD1 = "password"
const val PASSWORD2 = "password_2"

class IdGeneratorTest : StringSpec() {
    private val generator = IdGenerator(PASSWORD1)
    private val generator2 = IdGenerator(PASSWORD2)

    init {
        "same password, same input generates same ids" {
            forAll { a: String ->
                generator.generateIds(listOf(a)) == generator.generateIds(listOf(a))
            }
        }

        "different passwords generate different ids" {
            forNone { a: String ->
                generator.generateIds(listOf(a)) == generator2.generateIds(listOf(a))
            }
        }

        "repeated samples get one id" {
            forAll { a: String, b: String ->
                val samples = listOf(a, b, b, a)
                val ids = generator.generateIds(samples).values
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
                    generator.generateIds(listOf(a)) == generator.generateIds(listOf(b))
                } else {
                    generator.generateIds(listOf(a)) != generator.generateIds(listOf(b))
                }
            }
        }

        "ids are stable" {
            forAll { a: String, b: String, c: String ->
                val samples = listOf(a, b, c)
                val oldIds = generator.generateIds(samples).values
                val newIds = generator2.updateIds(PASSWORD1, samples, oldIds)
                samples.all {
                    (hmfId(it, generator, oldIds) != null) and (hmfId(it, generator, oldIds) == hmfId(it, generator2, newIds.values))
                }
            }
        }

        "ids for new patients are incremental" {
            forAll { a: String, b: String, c: String ->
                val samples = listOf(a, b)
                val newSamples = samples + c
                val oldIds = generator.generateIds(samples).values
                val newIds = generator2.updateIds(PASSWORD1, newSamples, oldIds)
                val previousMaxId = oldIds.maxBy { it.id }!!.id
                if (samples.contains(c)) {
                    if (c == a) {
                        hmfId(c, generator2, newIds.values)!! == 1
                    } else {
                        hmfId(c, generator2, newIds.values)!! == 2
                    }
                } else {
                    hmfId(c, generator2, newIds.values) == (previousMaxId + 1)
                }
            }
        }
    }

    private fun hmfId(sample: String, generator: IdGenerator, ids: Collection<HmfId>): Int? {
        val sampleHash = generator.generateIds(listOf(sample)).values.first().hash
        return ids.associateBy { it.hash }[sampleHash]?.id
    }
}
