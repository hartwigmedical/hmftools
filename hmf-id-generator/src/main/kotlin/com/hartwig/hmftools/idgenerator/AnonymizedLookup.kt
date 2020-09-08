package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.idgenerator.anonymizedIds.HmfSampleIdSimple

class AnonymizedLookup(password: String, private val sampleMap: Map<Hash, HmfSampleIdSimple>) {
    val generator = IdGenerator(password)

    companion object {
        operator fun invoke(password: String, samples: List<HmfSampleIdSimple>): AnonymizedLookup {
            return AnonymizedLookup(password, samples.map { Pair(it.hash, it) }.toMap())
        }
    }

    operator fun get(sample: String): HmfSampleIdSimple {
        return sampleMap[generator.hash(sample)]!!
    }
}