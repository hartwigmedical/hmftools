package com.hartwig.hmftools.idgenerator.anonymizedIds

import com.hartwig.hmftools.idgenerator.Hash
import com.hartwig.hmftools.idgenerator.IdGenerator

class AnonymizedMap(private val sampleMap: Map<String, Hash>, private val oldToNewMap: Map<Hash, Hash>) {

    companion object {
        fun create(oldPassword: String, newPassword: String, samples: List<String>): AnonymizedMap {
            val oldGenerator = IdGenerator(oldPassword)
            val newGenerator = IdGenerator(newPassword)

            val sampleMap = mutableMapOf<String, Hash>()
            val oldToNewMap = mutableMapOf<Hash, Hash>()
            for (sample in samples) {
                val oldHash = oldGenerator.hash(sample)
                val newHash = newGenerator.hash(sample)

                sampleMap[sample] = newHash
                oldToNewMap[oldHash] = newHash
            }

            return AnonymizedMap(sampleMap, oldToNewMap)
        }
    }

    fun oldHashes(): Set<Hash> {
        return oldToNewMap.keys
    }

    fun fromSample(sample: String): Hash {
        return sampleMap[sample]!!
    }

    fun fromOldHash(oldHash: Hash): Hash {
        return oldToNewMap[oldHash]!!
    }
}