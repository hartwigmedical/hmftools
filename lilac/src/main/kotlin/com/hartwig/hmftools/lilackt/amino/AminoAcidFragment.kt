package com.hartwig.hmftools.lilackt.amino

import com.hartwig.hmftools.lilackt.nuc.NucleotideFragment

class AminoAcidFragment(id: String, genes: Set<String>, nucleotideLoci: List<Int>, nucleotideQuality: List<Int>, nucleotides: List<String>,
                        private val aminoAcidLoci: List<Int>, private val aminoAcids: List<String>) : NucleotideFragment(id, genes, nucleotideLoci, nucleotideQuality, nucleotides) {


    fun containsAll(indices: Collection<Int>): Boolean {
        return indices.all { this.containsAminoAcid(it)}
    }

    fun containsAminoAcid(index: Int): Boolean {
        return aminoAcidLoci.contains(index)
    }

    fun aminoAcid(loci: Int): String {
        return aminoAcids[aminoAcidLoci.indexOf(loci)]
    }

    fun aminoAcids(vararg indices: Int): String {
        return indices.joinToString("") { aminoAcid(it) }
    }

    fun aminoAcidLoci(): List<Int> = aminoAcidLoci

    fun intersectAminoAcidLoci(otherAminoAcidLoci: Collection<Int>): AminoAcidFragment {
        val filteredIndexes = aminoAcidLoci
                .mapIndexed { index: Int, loci: Int -> Pair(index, loci) }
                .filter { (_, loci) -> loci in otherAminoAcidLoci }
                .map { (index, _) -> index }

        val filteredAminoAcidLoci = filteredIndexes.map { aminoAcidLoci[it] }
        val filteredAminoAcids = filteredIndexes.map { aminoAcids[it] }

        return AminoAcidFragment(id, genes, nucleotideLoci, nucleotideQuality, nucleotides, filteredAminoAcidLoci, filteredAminoAcids)
    }
}