package com.hartwig.hmftools.lilac.prot

import com.hartwig.hmftools.lilac.hla.HlaAllele

data class ProteinSequence(val contig: String, val proteins: String) {
    val allele: HlaAllele by lazy { HlaAllele(contig) }

    val fullProteinSequence: String by lazy { proteins.replace("*", "").replace(".", "") }
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
                val end = Math.min(proteins.length, boundary)
                result.add(proteins.substring(start, end).replace(".", "").replace("*", ""))
            }

            previousBoundary = boundary
        }

        return result
    }

    fun allKmers(length: Int): Set<String> {
        return fixedLengthRollingKmers(length, fullProteinSequence).toSet()
    }

    fun uniqueExonicKmers(length: Int, exonicBoundaries: List<Int>): Set<String> {
        return exonicProteins(exonicBoundaries)
                .flatMap { fixedLengthRollingKmers(length, it) }
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

    private fun variableLengthRollingKmers(minKmerSize: Int, sequence: String): List<String> {
        val result = mutableListOf<String>()
        for (i in 0..sequence.length - minKmerSize) {
            for (j in i + minKmerSize..sequence.length) {
                result.add(sequence.substring(i, j))
            }
        }

        return result
    }

}
