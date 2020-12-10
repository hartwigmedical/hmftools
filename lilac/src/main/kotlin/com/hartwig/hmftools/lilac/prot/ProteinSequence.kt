package com.hartwig.hmftools.lilac.prot

data class ProteinSequence(val contig: String, val proteins: String) {
    fun addProtein(more: String): ProteinSequence {
        return ProteinSequence(contig, proteins + more)
    }

    fun length(): Int {
        return removeSpecialCharacters().length
    }

    private fun removeSpecialCharacters(): String {
        return proteins.replace("*", "").replace(".", "")
    }

    fun exonicProteins(exonicBoundaries: List<Int>): List<String> {
        val result = mutableListOf<String>()
        var previousBoundary = -1
        for (i in exonicBoundaries.indices) {
            val boundary = exonicBoundaries[i]
            result.add(proteins.substring(previousBoundary + 1, boundary).replace(".", "").replace("*", ""))
            previousBoundary = boundary
        }

        return result
    }

    fun allKmers(length: Int): Set<String> {
        return rollingKmers(length, removeSpecialCharacters()).toSet()
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
