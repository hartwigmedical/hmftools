package com.hartwig.hmftools.lilac.prot

import com.hartwig.hmftools.lilac.hla.HlaAllele

data class ProteinSequence(val contig: String, val proteins: String) {
    val allele: HlaAllele by lazy { HlaAllele(contig) }

    private val fullProteinSequence: String by lazy { proteins.replace("*", "").replace(".", "") }
    val length: Int by lazy { fullProteinSequence.length }

    fun copyWithAdditionalProtein(more: String): ProteinSequence {
        return ProteinSequence(contig, proteins + more)
    }

    fun exonicProteins(exonicBoundaries: List<Int>): List<String> {
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
        return rollingKmers(length, fullProteinSequence).toSet()
    }

    fun exonicKmers(length: Int, exonicBoundaries: List<Int>): Set<String> {
        return exonicProteins(exonicBoundaries).flatMap { rollingKmers(length, it) }.toSet()
    }

    private fun rollingKmers(kmerSize: Int, sequence: String): List<String> {
        val result = mutableListOf<String>()
        for (i in 0..sequence.length - kmerSize) {
            result.add(sequence.substring(i, i + kmerSize))
        }


        return result
    }

}
