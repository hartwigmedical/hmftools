package com.hartwig.hmftools.lilac.prot

import com.hartwig.hmftools.lilac.hla.HlaAllele
import java.util.*

data class ProteinSequence(val contig: String, val proteins: String) {
    val allele: HlaAllele by lazy { HlaAllele(contig) }

    private val fullProteinSequence: String by lazy { proteins.removeSpecialCharacters() }
    val length: Int by lazy { fullProteinSequence.length }


    fun copyWithAdditionalProtein(more: String): ProteinSequence {
        return ProteinSequence(contig, proteins + more)
    }

    fun exonicProteins(exonicBoundaries: List<Int>): List<String> {
        if (exonicBoundaries.isEmpty()) {
            return listOf(fullProteinSequence)
        }

        val result = mutableListOf<String>()
        var previousBoundary = -1
        for (i in exonicBoundaries.indices) {
            val boundary = exonicBoundaries[i]
            val start = previousBoundary + 1
            if (start < proteins.length) {
                val end = proteins.length.coerceAtMost(boundary)
                result.add(proteins.substring(start, end).removeSpecialCharacters())
            }

            previousBoundary = boundary
        }

        if (previousBoundary < proteins.length - 1) {
            result.add(proteins.substring(previousBoundary + 1).removeSpecialCharacters())
        }

        return result
    }

    fun allKmers(length: Int): Set<String> {
        return fixedLengthRollingKmers(length, fullProteinSequence).toSet()
    }

    fun uniqueExonicKmers(kmerLength: Int, exonicBoundaries: List<Int>): Set<String> {
        return uniqueExonicKmers(kmerLength, kmerLength, exonicBoundaries)
    }

    fun uniqueExonicKmers(minKmerLength: Int, maxKmerLength: Int, exonicBoundaries: List<Int>): Set<String> {
        return exonicProteins(exonicBoundaries)
                .flatMap { variableLengthRollingKmers(minKmerLength, maxKmerLength, it) }
                .groupBy { it }
                .filter { it.value.size == 1 }
                .keys
    }

    private fun fixedLengthRollingKmers(kmerSize: Int, sequence: String): List<String> {
        val result = mutableListOf<String>()
        for (i in 0..sequence.length - kmerSize) {
            result.add(sequence.substring(i, i + kmerSize))
        }

        return result
    }

    private fun variableLengthRollingKmers(minKmerLength: Int, maxKmerLength: Int, sequence: String): List<String> {
        if (sequence.length < minKmerLength) {
            return Collections.emptyList()
        }

        return fixedLengthRollingKmers(maxKmerLength.coerceAtMost(sequence.length), sequence);
    }

    private fun String.removeSpecialCharacters(): String {
        return this.replace("*", "").replace(".", "");
    }

}
