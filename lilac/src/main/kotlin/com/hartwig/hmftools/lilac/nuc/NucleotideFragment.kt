package com.hartwig.hmftools.lilac.nuc

import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.lilac.amino.AminoAcidFragment
import com.hartwig.hmftools.lilac.read.SAMRecordRead

open class NucleotideFragment(
        val id: String,
        val genes: Set<String>,
        protected val nucleotideLoci: List<Int>,
        protected val nucleotideQuality: List<Int>,
        protected val nucleotides: List<Char>) {

    companion object {

        fun fromReads(minBaseQual: Int, reads: List<SAMRecordRead>): List<NucleotideFragment> {
            return reads.groupBy { it.samRecord.readName }.map { fromReadPairs(minBaseQual, it.value) }
        }

        private fun fromReadPairs(minBaseQual: Int, reads: List<SAMRecordRead>): NucleotideFragment {
            fun nucleotide(index: Int): Char {
                for (read in reads) {
                    if (read.nucleotideIndices().contains(index)) {
                        return read.nucleotide(index)
                    }
                }
                throw IllegalArgumentException("Fragment does not contain nucleotide at location $index")
            }

            fun quality(index: Int): Int {
                for (read in reads) {
                    if (read.nucleotideIndices().contains(index)) {
                        return read.quality(index)
                    }
                }
                throw IllegalArgumentException("Fragment does not contain nucleotide at location $index")
            }

            val nucleotideIndices = reads
                    .flatMap { it.nucleotideIndices() }
                    .distinct()
                    .sorted()

            val nucleotides = nucleotideIndices.map { nucleotide(it) }
            val nucleotideQuality = nucleotideIndices.map { quality(it) }

            val id = reads[0].samRecord.readName + " -> " + reads[0].samRecord.alignmentStart
            val genes = mutableSetOf<String>()
            genes.add(reads[0].gene)
            if (reads.size > 1) {
                genes.add(reads[1].gene)
            }

            return NucleotideFragment(id, genes, nucleotideIndices, nucleotideQuality, nucleotides)
        }
    }

    fun isEmpty(): Boolean {
        return nucleotideLoci.isEmpty()
    }

    fun isNotEmpty(): Boolean {
        return !isEmpty()
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

    fun nucleotide(loci: Int): Char {
        return nucleotides[nucleotideLoci.indexOf(loci)]
    }

    fun nucleotides(): List<Char> = nucleotides

    fun nucleotideLoci(): List<Int> = nucleotideLoci

    fun nucleotideQuality(): List<Int> = nucleotideQuality

    fun qualityFilter(minBaseQual: Int): NucleotideFragment {
        val qualityFilteredIndexes = nucleotideQuality
                .mapIndexed { index: Int, quality: Int -> Pair(index, quality) }
                .filter { (_, quality) -> quality >= minBaseQual }
                .map { (index, _) -> index }

        val qualityFilteredNucleotideLoci = qualityFilteredIndexes.map { nucleotideLoci[it] }
        val qualityFilteredNucleotides = qualityFilteredIndexes.map { nucleotides[it] }
        val qualityFilteredNucleotideQuality = qualityFilteredIndexes.map { nucleotideQuality[it] }

        return NucleotideFragment(id, genes, qualityFilteredNucleotideLoci, qualityFilteredNucleotideQuality, qualityFilteredNucleotides)
    }

    fun toAminoAcidFragment(): AminoAcidFragment {
        fun aminoAcid(index: Int): Char {
            val first = nucleotide(index * 3)
            val second = nucleotide(index * 3 + 1)
            val third = nucleotide(index * 3 + 2)
            return Codons.aminoAcid(first.toString() + second + third)
        }

        val aminoAcidIndices = nucleotideLoci
                .filter { it % 3 == 0 }
                .filter { nucleotideLoci.contains(it + 1) && nucleotideLoci.contains(it + 2) }
                .map { it / 3 }

        val aminoAcids = aminoAcidIndices.map { aminoAcid(it) }

        return AminoAcidFragment(id, genes, nucleotideLoci, nucleotideQuality, nucleotides, aminoAcidIndices, aminoAcids)
    }

    fun enrich(loci: Int, nucleotide: Char, quality: Int): NucleotideFragment {
        assert(!containsNucleotide(loci))
        return NucleotideFragment(id, genes, nucleotideLoci + loci, nucleotideQuality + quality, nucleotides + nucleotide)
    }

}