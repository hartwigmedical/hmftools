package com.hartwig.hmftools.lilac.prot

import com.hartwig.hmftools.lilac.ext.rollingKmers
import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.phase.PhasedEvidence2

data class ProteinSequence(val contig: String, val protein: String) {
    val allele: HlaAllele by lazy { HlaAllele(contig) }

    val fullProteinSequence: String by lazy { protein.removeIndels() }
    val length: Int by lazy { fullProteinSequence.length }


    fun copyWithAdditionalProtein(more: String): ProteinSequence {
        return ProteinSequence(contig, protein + more)
    }


    fun consistentWith(evidence: PhasedEvidence2): Boolean {
        return evidence.evidence.keys.any { this.consistentWith(evidence.aminoAcidIndices, it.toCharArray()) }
    }

    fun consistentWith(aminoAcidIndices: IntArray, aminoAcids: CharArray): Boolean {
        for (i in aminoAcidIndices.indices) {
            val index = aminoAcidIndices[i]
            val aminoAcid = aminoAcids[i]

            if (this.length > index && this.fullProteinSequence[index] != '*' && this.fullProteinSequence[index] != aminoAcid) {
                return false
            }
        }

        return true
    }

    fun inflate(template: String): ProteinSequence {
        val joiner = StringBuilder()
        for (i in protein.indices) {
            if (protein[i] == '-') {
                joiner.append(template[i])
            } else {
                joiner.append(protein[i])
            }
        }
        return ProteinSequence(contig, joiner.toString())
    }

    fun deflate(template: String): ProteinSequence {
        val joiner = StringBuilder()
        for (i in protein.indices) {
            when {
                protein[i] == '.' -> joiner.append('.')
                protein[i] == template[i] -> joiner.append('-')
                else -> joiner.append(protein[i])
            }
        }
        return ProteinSequence(contig, joiner.toString())
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
            if (start < protein.length) {
                val end = protein.length.coerceAtMost(boundary)
                result.add(protein.substring(start, end).removeSpecialCharacters())
            }

            previousBoundary = boundary
        }

        if (previousBoundary < protein.length - 1) {
            result.add(protein.substring(previousBoundary + 1).removeSpecialCharacters())
        }

        return result.filter { it.isNotEmpty() }
    }

    fun allKmers(length: Int): Set<String> {
        return fullProteinSequence.rollingKmers(length).toSet()
    }

    fun exonicKmers(kmerLength: Int, exonicBoundaries: List<Int>): Set<String> {
        return exonicKmers(kmerLength, kmerLength, exonicBoundaries)
    }

    fun exonicKmers(minKmerLength: Int, maxKmerLength: Int, exonicBoundaries: List<Int>): Set<String> {
        return exonicProteins(exonicBoundaries)
                .flatMap { it.rollingKmers(minKmerLength, maxKmerLength) }
                .toSet()
//                .groupBy { it }
//                .filter { it.value.size == 1 }
//                .keys
    }


    private fun String.removeIndels(): String {
        return this.replace(".", "");
    }

    private fun String.removeSpecialCharacters(): String {
        return this.replace("*", "").replace(".", "");
    }

}
