package com.hartwig.hmftools.lilac.read

import com.hartwig.hmftools.lilac.nuc.NucleotideFragment

class Fragment(alignedGene: String, id: String, nucleotideLoci: List<Int>, nucleotides: List<Char>, private val aminoAcidLoci: List<Int>, private val aminoAcids: List<Char>) :
        NucleotideFragment(alignedGene, id, nucleotideLoci, nucleotides) {

    fun containsAminoAcid(index: Int): Boolean {
        return aminoAcidLoci.contains(index)
    }

    fun aminoAcid(loci: Int): Char {
        return aminoAcids[aminoAcidLoci.indexOf(loci)]
    }

    fun aminoAcidIndices(): List<Int> = aminoAcidLoci
}