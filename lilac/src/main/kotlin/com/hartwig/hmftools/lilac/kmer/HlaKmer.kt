package com.hartwig.hmftools.lilac.kmer

import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.prot.ProteinSequence

class HlaKmer(a: Map<ProteinSequence, Set<String>>, b: Map<ProteinSequence, Set<String>>, c: Map<ProteinSequence, Set<String>>) {
    private val kmerCount = mutableMapOf<String, Int>()
    private val sequences = mutableMapOf<ProteinSequence, Set<String>>()

    init {
        sequences.putAll(a)
        sequences.putAll(b)
        sequences.putAll(c)
        sequences.values.flatten().forEach { kmer -> kmerCount.compute(kmer) { _, old -> (old ?: 0) + 1 } }
    }

    val uniqueKmers = kmerCount.filter { it.value == 1 }.map { it.key }.toSet()


    fun proteinSequence(uniqueKmer: String): ProteinSequence {
        assert(uniqueKmers.contains(uniqueKmer))

        val aLookup = sequences.filter { it.value.contains(uniqueKmer) }.map { it.key }
        if (aLookup.isNotEmpty()) {
            return aLookup[0]
        }

        throw IllegalArgumentException()
    }


    fun uniqueKmers(): Set<String> {
        return uniqueKmers;
    }

    fun kmers(): Set<String> {
        return kmerCount.keys
    }

    fun kmers(hlaAlleles: Set<HlaAllele>): Set<String> {
        return sequences(hlaAlleles).map { it.value }.flatten().toSet()
    }

    fun sequences(hlaAlleles: Set<HlaAllele>): Map<ProteinSequence, Set<String>> {
        return sequences.filter { hlaAlleles.contains(it.key.allele) }
    }

    fun sequences(): Map<ProteinSequence, Set<String>> {
        return sequences
    }


}