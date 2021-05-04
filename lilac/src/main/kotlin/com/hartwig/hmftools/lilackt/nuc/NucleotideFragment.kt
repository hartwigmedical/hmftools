package com.hartwig.hmftools.lilackt.nuc

import com.hartwig.hmftools.common.codon.Codons
import com.hartwig.hmftools.lilackt.amino.AminoAcidFragment

open class NucleotideFragment(
        val id: String,
        val genes: Set<String>,
        protected val nucleotideLoci: List<Int>,
        protected val nucleotideQuality: List<Int>,
        protected val nucleotides: List<String>) {

    companion object {

        fun merge(o1: NucleotideFragment, o2: NucleotideFragment): NucleotideFragment {
            require(o1.id == o2.id)
            val genes = o1.genes union o2.genes
            return NucleotideFragment(o1.id, genes, o1.nucleotideLoci + o2.nucleotideLoci, o1.nucleotideQuality + o2.nucleotideQuality, o1.nucleotides + o2.nucleotides)
        }
    }


    fun isEmpty(): Boolean {
        return nucleotideLoci.isEmpty()
    }

    fun isNotEmpty(): Boolean {
        return !isEmpty()
    }

    fun containsIndel(): Boolean {
        return nucleotides.any { it == "." || it.length > 1 }
    }

    fun containsNucleotide(index: Int): Boolean {
        return nucleotideLoci.contains(index)
    }

    fun containsAllNucleotides(vararg indices: Int): Boolean {
        return indices.all { containsNucleotide(it) }
    }

    fun nucleotides(vararg indices: Int): String {
        return indices.joinToString("") { nucleotide(it) }
    }

    fun nucleotide(loci: Int): String {
        return nucleotides[nucleotideLoci.indexOf(loci)]
    }

    fun nucleotides(): List<String> = nucleotides

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
        fun aminoAcid(index: Int): String {
            val first = nucleotide(index * 3)
            val second = nucleotide(index * 3 + 1)
            val third = nucleotide(index * 3 + 2)
            if (first == "." && second == "." && third == ".") {
                return "."
            }

            return Codons.aminoAcids(first + second + third)
        }

        val aminoAcidIndices = nucleotideLoci
                .filter { it % 3 == 0 }
                .filter { nucleotideLoci.contains(it + 1) && nucleotideLoci.contains(it + 2) }
                .map { it / 3 }

        val aminoAcids = aminoAcidIndices.map { aminoAcid(it) }.map { it }

        return AminoAcidFragment(id, genes, nucleotideLoci, nucleotideQuality, nucleotides, aminoAcidIndices, aminoAcids)
    }

    fun enrich(loci: Int, nucleotide: String, quality: Int): NucleotideFragment {
        assert(!containsNucleotide(loci))
        return NucleotideFragment(id, genes, nucleotideLoci + loci, nucleotideQuality + quality, nucleotides + nucleotide)
    }

}