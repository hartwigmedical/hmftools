package com.hartwig.hmftools.idgenerator

import org.bouncycastle.jcajce.provider.digest.SHA3
import org.bouncycastle.util.encoders.Hex

class IdGenerator(private val password: String) {
    fun generateIds(patientIds: List<String>): Set<HmfId> {
        return patientIds.distinct().mapIndexed { index, patientId -> HmfId(hash(patientId), index + 1) }.toSet()
    }

    fun updateIds(oldPassword: String, patientIds: List<String>, oldIds: Set<HmfId>): Set<HmfId> {
        val indexedOldIds = oldIds.associateBy { it.hash }
        val oldGenerator = IdGenerator(oldPassword)
        val oldNewHashPairs = patientIds.distinct().map { Pair(oldGenerator.hash(it), hash(it)) }
        assertIncludesOldPatients(indexedOldIds, oldNewHashPairs)
        return updateIds(oldNewHashPairs, indexedOldIds)
    }

    private fun updateIds(hashPairs: List<Pair<String, String>>, indexedOldIds: Map<String, HmfId>): Set<HmfId> {
        val highestId = indexedOldIds.maxBy { it.value.id }?.value?.id ?: 0
        val (idsToUpdate, idsToGenerate) = hashPairs.partition { indexedOldIds.containsKey(it.first) }
        val updatedIds = idsToUpdate.map { HmfId(it.second, indexedOldIds[it.first]!!.id) }
        val newIds = idsToGenerate.mapIndexed { index, (_, newHash) -> HmfId(newHash, (highestId + 1 + index)) }
        return (updatedIds + newIds).toSet()
    }

    //MIVO: check that all old patients are included in the new list
    private fun assertIncludesOldPatients(oldIds: Map<String, HmfId>, hashPairs: List<Pair<String, String>>) {
        val hashes = hashPairs.map { it.first }.toSet()
        if (!oldIds.keys.all { hashes.contains(it) }) {
            throw IllegalArgumentException("New patient list does not include all previous patients")
        }
    }

    private fun hash(input: String): String {
        val sha3 = SHA3.Digest256()
        sha3.update((input + password).toByteArray())
        return Hex.toHexString(sha3.digest())
    }
}
