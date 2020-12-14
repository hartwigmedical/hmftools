package com.hartwig.hmftools.lilac.kmer

import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.common.utils.SuffixTree
import htsjdk.samtools.SAMRecord
import java.util.concurrent.ConcurrentHashMap
import java.util.function.Consumer

class BamKmer(private val codonKmers: Set<String>) : Consumer<SAMRecord> {
    private val map = ConcurrentHashMap<String, Int>()

    fun kmerCount(): Map<String, Int> {
        return map
    }

    override fun accept(t: SAMRecord) {
        for ((kmer, count) in kmerCountDna(t)) {
            map.compute(kmer) { _: String, y: Int? -> (y ?: 0) + count }
        }
    }

    private fun kmerCountDna(record: SAMRecord): Map<String, Int> {
        return kmerCountDna(record.readString)
    }

    private fun reverseComplement(dna: String): String {
        val builder = StringBuilder()
        for (i in dna.indices.reversed()) {
            when (dna[i]) {
                'G' -> builder.append('C')
                'C' -> builder.append('G')
                'A' -> builder.append('T')
                'T' -> builder.append('A')
                else -> builder.append(dna[i])

            }
        }


        return builder.toString()
    }

    private fun kmerCountDna(dna: String): Map<String, Int> {
        val dnaReversed = reverseComplement(dna)

        val codonString0 = Codons.asCodonString(dna)
        val codonString1 = Codons.asCodonString(dna.substring(1))
        val codonString2 = Codons.asCodonString(dna.substring(2))
        val codonString3 = Codons.asCodonString(dnaReversed)
        val codonString4 = Codons.asCodonString(dnaReversed.substring(1))
        val codonString5 = Codons.asCodonString(dnaReversed.substring(2))

        val searchTrees = listOf(
                SuffixTree(codonString0),
                SuffixTree(codonString1),
                SuffixTree(codonString2),
                SuffixTree(codonString3),
                SuffixTree(codonString4),
                SuffixTree(codonString5))

        val result = mutableMapOf<String, Int>()
        for (kmer in codonKmers) {
            for (tree in searchTrees) {
                if (tree.contains(kmer)) {
                    result.compute(kmer) { _: String, y: Int? -> (y ?: 0) + 1 }
                }
            }
        }

        return result
    }
}