package com.hartwig.hmftools.lilac.read

import com.hartwig.hmftools.common.codon.Codons

data class Fragment(val reads: List<Read>) {
    private val nucleotideIndices = reads
            .flatMap { it.nucleotideIndices().toList() }
            .distinct()
            .sorted()

    private val aminoAcidIndices = nucleotideIndices
            .filter { it % 3 == 0 }
            .filter { nucleotideIndices.contains(it + 1) && nucleotideIndices.contains(it + 2) }
            .map { it / 3 }

    fun containsAminoAcid(index: Int): Boolean {
        return aminoAcidIndices.contains(index)
    }

    fun aminoAcid(index: Int, minQual: Int): Char {
        val first = nucleotide(index * 3, minQual)
        if (first == '.') {
            return '.'
        }
        val second = nucleotide(index * 3 + 1, minQual)
        if (second == '.') {
            return '.'
        }
        val third = nucleotide(index * 3 + 2, minQual)
        if (third == '.') {
            return '.'
        }

        return Codons.aminoAcid(first.toString() + second + third)
    }

    fun nucleotide(index: Int, minQual: Int): Char {
        for (read in reads) {
            if (read.nucleotideIndices().contains(index)) {
                return read.nucleotide(index, minQual)
            }
        }

        throw IllegalArgumentException("Fragment does not contain nucleotide at location $index")
    }

    fun aminoAcidIndices(): List<Int> = aminoAcidIndices

    fun nucleotideIndices(): List<Int> = nucleotideIndices

}