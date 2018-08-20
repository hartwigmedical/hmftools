package com.hartwig.hmftools.idgenerator

import org.bouncycastle.jcajce.provider.digest.SHA3
import org.bouncycastle.util.encoders.Hex
import kotlin.math.sqrt

typealias CpctId = String
typealias OldHash = String
typealias NewHash = String

class IdGenerator(private val password: String) {
    private data class HashTriple(val cpctId: CpctId, val oldHash: OldHash, val newHash: NewHash)

    fun generateIds(patientIds: List<CpctId>): Map<CpctId, HmfId> {
        return patientIds.distinct().mapIndexed { index, patientId ->
            Pair(patientId, HmfId(hash(patientId), index + 1))
        }.toMap()
    }

    fun updateIds(oldPassword: String, patientIds: List<CpctId>, oldIds: Collection<HmfId>): Map<CpctId, HmfId> {
        val idsPerOldHash = oldIds.distinct().associateBy { it.hash }
        val oldGenerator = IdGenerator(oldPassword)
        val hashTriples = patientIds.distinct().map { HashTriple(it, oldGenerator.hash(it), hash(it)) }
        checkOldPasswordMatched(idsPerOldHash, hashTriples)
        return updateIds(hashTriples, idsPerOldHash)
    }

    fun hash(input: String): String {
        val sha3 = SHA3.Digest256()
        sha3.update((input + password).toByteArray())
        return Hex.toHexString(sha3.digest())
    }

    private fun updateIds(hashTriples: List<HashTriple>, idsPerOldHash: Map<OldHash, HmfId>): Map<CpctId, HmfId> {
        val highestId = idsPerOldHash.maxBy { it.value.id }?.value?.id ?: 0
        val (idsToUpdate, idsToGenerate) = hashTriples.partition { idsPerOldHash.containsKey(it.oldHash) }
        val updatedIds = idsToUpdate.map { updateHmfId(it, getOldId(idsPerOldHash, it.oldHash)) }
        val newIds = idsToGenerate.mapIndexed { index, triple -> updateHmfId(triple, (highestId + 1 + index)) }
        return (updatedIds + newIds).toMap()
    }

    //MIVO: check that at least some of the old patients are included in the new list
    private fun checkOldPasswordMatched(oldIds: Map<OldHash, HmfId>, hashTriples: List<HashTriple>) {
        val oldHashMatches = numberOfOldHashMatches(oldIds, hashTriples)
        if (oldHashMatches <= sqrt(oldIds.size.toDouble())) {
            error("Could only reproduce $oldHashMatches of the existing hashes using the provided $OLD_PASSWORD and $PATIENT_IDS_FILE. " +
                          "Either too many patients were removed from the $PATIENT_IDS_FILE file or $OLD_PASSWORD was wrong.")
        }
    }

    private fun numberOfOldHashMatches(oldIds: Map<OldHash, HmfId>, hashTriples: List<HashTriple>): Int {
        val oldPasswordHashes = hashTriples.map { it.oldHash }.toSet()
        return oldPasswordHashes.filter { oldIds.containsKey(it) }.size
    }

    private fun getOldId(idsPerOldHash: Map<OldHash, HmfId>, oldHash: OldHash): Int {
        return idsPerOldHash[oldHash]!!.id
    }

    private fun updateHmfId(hashTriple: HashTriple, newId: Int): Pair<CpctId, HmfId> {
        return Pair(hashTriple.cpctId, HmfId(hashTriple.newHash, newId))
    }
}
