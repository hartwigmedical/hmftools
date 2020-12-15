package com.hartwig.hmftools.lilac.kmer

import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.prot.ProteinSequence
import org.apache.logging.log4j.LogManager

class HlaKmer private constructor(val sequences: Map<ProteinSequence, Set<String>>) {

    companion object {
        val logger = LogManager.getLogger(this::class.java)

        operator fun invoke(collection: Collection<Map<ProteinSequence, Set<String>>>): HlaKmer {
            var excluded = 0
            val combined = mutableMapOf<ProteinSequence, Set<String>>()
            for (map in collection) {
                for ((proteinSequence, kmers) in map.entries) {
                    if (combined.containsValue(kmers)) {
                        logger.debug("Excluding indistinguishable hla type ${proteinSequence.allele}")

                        excluded++
                    } else {
                        combined[proteinSequence] = kmers
                    }
                }
            }

            if (excluded > 0) {
                logger.warn("Excluded $excluded indistinguishable hla types")
            }

            return HlaKmer(combined)
        }
    }

    val kmers = sequences.values.flatten().toSet()
    val uniqueKmers = sequences.values.flatten().groupBy { it }.filter { it.value.size == 1 }.keys

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
        return kmers
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