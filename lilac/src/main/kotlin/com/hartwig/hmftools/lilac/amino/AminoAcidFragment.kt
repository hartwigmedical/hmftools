package com.hartwig.hmftools.lilac.amino

import com.hartwig.hmftools.lilac.nuc.NucleotideFragment

class AminoAcidFragment(id: String, genes: Set<String>, nucleotideLoci: List<Int>, nucleotideQuality: List<Int>, nucleotides: List<Char>,
                        private val aminoAcidLoci: List<Int>, private val aminoAcids: List<Char>) : NucleotideFragment(id, genes, nucleotideLoci, nucleotideQuality, nucleotides) {

    fun containsAminoAcid(index: Int): Boolean {
        return aminoAcidLoci.contains(index)
    }

    fun aminoAcid(loci: Int): Char {
        return aminoAcids[aminoAcidLoci.indexOf(loci)]
    }

    fun aminoAcidIndices(): List<Int> = aminoAcidLoci

    fun intersectAminoAcidLoci(otherAminoAcidLoci: Collection<Int>): AminoAcidFragment {
        val filteredIndexes = aminoAcidLoci
                .mapIndexed { index: Int, loci: Int -> Pair(index, loci) }
                .filter { (_, loci) -> loci in otherAminoAcidLoci }
                .map { (index, _) -> index }

        val filteredAminoAcidLoci = filteredIndexes.map { aminoAcidLoci[it] }
        val filteredAminoAcids = filteredIndexes.map { aminoAcids[it] }

        return AminoAcidFragment(id, genes, nucleotideLoci, nucleotideQuality, nucleotides, filteredAminoAcidLoci, filteredAminoAcids)
    }

    fun qualityFilterNucleotides(minBaseQuality: Int): AminoAcidFragment {
        val qualityFilteredNucleoties = this.qualityFilter(minBaseQuality)
        return AminoAcidFragment(id,
                genes,
                qualityFilteredNucleoties.nucleotideLoci(),
                qualityFilteredNucleoties.nucleotideQuality(),
                qualityFilteredNucleoties.nucleotides(),
                aminoAcidLoci,
                aminoAcids)
    }

}