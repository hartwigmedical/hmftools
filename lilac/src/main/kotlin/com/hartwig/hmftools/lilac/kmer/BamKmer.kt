package com.hartwig.hmftools.lilac.kmer

import com.hartwig.hmftools.common.utils.SuffixTree
import com.hartwig.hmftools.lilac.dna.aminoAcids
import com.hartwig.hmftools.lilac.dna.dnaReverseComplement
import htsjdk.samtools.SAMRecord
import java.util.concurrent.ConcurrentHashMap
import java.util.function.Consumer

class BamKmer(private val codonKmers: Set<String>) : Consumer<SAMRecord> {
    private val map = ConcurrentHashMap<String, Int>()

    fun kmerCount(): Map<String, Int> {
        return map
    }

    override fun accept(t: SAMRecord) {

        val counts = kmerCountDna(t)
        for ((kmer, count) in counts) {
            map.compute(kmer) { _: String, y: Int? -> (y ?: 0) + count }
        }
    }

    private fun kmerCountDna(record: SAMRecord): Map<String, Int> {
        return kmerCountDna(record.readString)
    }


    private fun kmerCountDna(dna: String): Map<String, Int> {
        val dnaReversed = dna.dnaReverseComplement()

        val codonString0 = dna.aminoAcids()
        val codonString1 = dna.substring(1).aminoAcids()
        val codonString2 = dna.substring(2).aminoAcids()
        val codonString3 = dnaReversed.aminoAcids()
        val codonString4 = dnaReversed.substring(1).aminoAcids()
        val codonString5 = dnaReversed.substring(2).aminoAcids()

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