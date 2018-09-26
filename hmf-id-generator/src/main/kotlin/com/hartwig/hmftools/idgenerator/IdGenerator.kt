package com.hartwig.hmftools.idgenerator

import com.hartwig.hmftools.idgenerator.anonymizedIds.HashId
import org.bouncycastle.jcajce.provider.digest.SHA3
import org.bouncycastle.util.encoders.Hex


class IdGenerator(private val password: String, private val defaultSeed: Int = 0) {

    fun generate(input: Collection<String>, seed: Int = 0): Collection<HashId> {
        return input.distinct().mapIndexed { index, string -> HashId(hash(string), index + 1 + seed) }
    }

    fun update(newPassword: String, input: Collection<String>, existingIds: Collection<HashId>): Collection<HashId> {
        val updatedIds = updateIdList(input, existingIds)
        return changePassword(newPassword, input, updatedIds)
    }

    //MIVO: change current password
    fun changePassword(newPassword: String, input: Collection<String>, existingIds: Collection<HashId>): Collection<HashId> {
        checkPasswordMatchesForAllInputs(input, existingIds)
        checkAllPlaintextsExist(input, existingIds)
        val existingIdMap = existingIds.associateBy { it.hash }
        val newGenerator = IdGenerator(newPassword)
        val updatedIds = input.distinct().map { string -> existingIdMap[hash(string)]!!.copy(hash = newGenerator.hash(string)) }
        return removedIds(input, existingIds) + updatedIds
    }

    fun updateIdList(input: Collection<String>, existingIds: Collection<HashId>): Collection<HashId> {
        checkPasswordMatchesForSomeInputs(input, existingIds)
        val existingHashes = existingIds.map { it.hash }.toSet()
        val idsToGenerate = input.filterNot { existingHashes.contains(hash(it)) }.toSet()
        return existingIds + generate(idsToGenerate, highestId(existingIds))
    }

    fun hash(plaintext: String): Hash {
        val sha3 = SHA3.Digest256()
        sha3.update((plaintext + password).toByteArray())
        return Hash(Hex.toHexString(sha3.digest()))
    }

    private fun highestId(existingIds: Collection<HashId>): Int = existingIds.maxBy { it.id }?.id ?: defaultSeed

    private fun removedIds(input: Collection<String>, existingIds: Collection<HashId>): Collection<HashId> {
        val inputHashes = input.map { hash(it) }.toSet()
        return existingIds.filterNot { inputHashes.contains(it.hash) }
    }

    //MIVO: check that the old password is correct for all strings
    private fun checkPasswordMatchesForAllInputs(input: Collection<String>, existingIds: Collection<HashId>) {
        val hashes = input.map { hash(it) }.toSet()
        val existingHashes = existingIds.map { it.hash }.toSet()
        if (hashes.isNotEmpty() && hashes.any { !existingHashes.contains(it) }) {
            error("$NEW_PASSWORD does not match")
        }
    }

    //MIVO: check that plaintext collection matches the existing id collection
    private fun checkAllPlaintextsExist(strings: Collection<String>, existingIds: Collection<HashId>) {
        val hashes = strings.map { hash(it) }.toSet()
        if (hashes.any { !existingIds.map { it.hash }.contains(it) }) {
            println(hashes)
            println(existingIds.map { it.hash }.toSet())
            error("Current list of strings to be hashed does not match existing list.")
        }
    }

    //MIVO: check that the old password matches for at least some of the strings
    private fun checkPasswordMatchesForSomeInputs(input: Collection<String>, existingIds: Collection<HashId>) {
        val hashes = input.map { hash(it) }.toSet()
        val existingHashes = existingIds.map { it.hash }.toSet()
        val numMatches = hashes.count { existingHashes.contains(it) }
        if (hashes.isNotEmpty() && existingHashes.isNotEmpty() && numMatches < 1) {
            error("Could not reproduce any of the existing hashes using the provided $NEW_PASSWORD and $SAMPLE_IDS_FILE. " +
                          "Most likely $NEW_PASSWORD was wrong. ")
        }
    }
}
