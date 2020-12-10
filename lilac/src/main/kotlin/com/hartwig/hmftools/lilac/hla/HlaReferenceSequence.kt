package com.hartwig.hmftools.lilac.hla

import htsjdk.samtools.reference.FastaSequenceFile
import htsjdk.samtools.reference.ReferenceSequence
import org.apache.logging.log4j.LogManager
import java.io.File

data class HlaReferenceSequence(val allele: HlaAllele, val length: Int, val name: String, val sequence: String) {

    companion object {

        val logger = LogManager.getLogger(this::class.java)

        operator fun invoke(file: File): List<HlaReferenceSequence> {
            val map = mutableMapOf<String, HlaReferenceSequence>()

            FastaSequenceFile(file, false).use {

                var nextSequence = it.nextSequence();
                while (nextSequence != null) {
                    val newSequence = HlaReferenceSequence(nextSequence)
                    val fourDigitName = newSequence.allele.fourDigitName()
                    if (!map.containsKey(fourDigitName)) {
                        map[fourDigitName] = newSequence
                    } else {
                        val existing = map[fourDigitName]!!
                        if (existing.length < newSequence.length) {
                            logger.info("Replacing ${existing.allele} with ${newSequence.allele}")
                            map[fourDigitName] = newSequence
                        }
                    }

                    nextSequence = it.nextSequence();
                }
            }

            return map.values.sortedBy { it.allele.fourDigitName() }
        }

        private operator fun invoke(referenceSequence: ReferenceSequence): HlaReferenceSequence {
            val name = referenceSequence.name
            val nameSplit = name.split(" ")
            val length = nameSplit[2].toInt()
            val contig = nameSplit[1]
            val allele = HlaAllele(contig)
            return HlaReferenceSequence(allele, length, name, referenceSequence.baseString)
        }
    }

    fun rollingKmers(kmerSize: Int): Set<String> {
        val result = mutableSetOf<String>()
        for (i in 0..sequence.length - kmerSize) {
            result.add(sequence.substring(i, i + kmerSize))
        }


        return result
    }

    fun containsKmer(kmer: String): Boolean {
        return sequence.contains(kmer)
    }

}