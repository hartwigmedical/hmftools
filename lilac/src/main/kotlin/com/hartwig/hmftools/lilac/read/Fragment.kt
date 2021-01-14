package com.hartwig.hmftools.lilac.read

import com.hartwig.hmftools.common.codon.Codons

class Fragment(val id: String, private val nucleotideLoci: List<Int>, private val nucleotides: List<Char>, private val aminoAcidLoci: List<Int>, private val aminoAcids: List<Char>) {

    companion object {

        fun fromReads(minBaseQual: Int, reads: List<SAMRecordRead>): List<Fragment> {
            return reads.groupBy { it.samRecord.readName }.map { fromReadPairs(minBaseQual, it.value) }
        }

        private fun fromReadPairs(minBaseQual: Int, reads: List<SAMRecordRead>): Fragment {
            fun nucleotide(index: Int): Char {
                for (read in reads) {
                    if (read.nucleotideIndices().contains(index)) {
                        return read.nucleotide(index)
                    }
                }
                throw IllegalArgumentException("Fragment does not contain nucleotide at location $index")
            }

            fun aminoAcid(index: Int): Char {
                val first = nucleotide(index * 3)
                val second = nucleotide(index * 3 + 1)
                val third = nucleotide(index * 3 + 2)
                return Codons.aminoAcid(first.toString() + second + third)
            }

            val nucleotideIndices = reads
                    .flatMap { it.nucleotideIndices(minBaseQual).toList() }
                    .distinct()
                    .sorted()

            val aminoAcidIndices = nucleotideIndices
                    .filter { it % 3 == 0 }
                    .filter { nucleotideIndices.contains(it + 1) && nucleotideIndices.contains(it + 2) }
                    .map { it / 3 }

            val nucleotides = nucleotideIndices.map { nucleotide(it) }
            val aminoAcids = aminoAcidIndices.map { aminoAcid(it) }
            val id = reads[0].samRecord.readName + " -> " + reads[0].samRecord.alignmentStart

            return Fragment(id, nucleotideIndices, nucleotides, aminoAcidIndices, aminoAcids)
        }

    }

    fun containsAminoAcid(index: Int): Boolean {
        return aminoAcidLoci.contains(index)
    }

    fun containsNucleotide(index: Int): Boolean {
        return nucleotideLoci.contains(index)
    }

    fun containsAllNucleotides(vararg indices: Int): Boolean {
        return indices.all { containsNucleotide(it) }
    }

    fun nucleotides(vararg indices: Int): String {
        return indices.map { nucleotide(it) }.joinToString("")
    }

    fun aminoAcid(loci: Int): Char {
        return aminoAcids[aminoAcidLoci.indexOf(loci)]
    }

    fun nucleotide(loci: Int): Char {
        return nucleotides[nucleotideLoci.indexOf(loci)]
    }

    fun aminoAcidIndices(): List<Int> = aminoAcidLoci

    fun nucleotideIndices(): List<Int> = nucleotideLoci

}