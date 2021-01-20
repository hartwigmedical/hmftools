package com.hartwig.hmftools.lilac.read

import com.hartwig.hmftools.lilac.nuc.NucleotideFragment

class Fragment(id: String, nucleotideLoci: List<Int>, nucleotideQuality: List<Int>, nucleotides: List<Char>, genes: Set<String>, private val aminoAcidLoci: List<Int>, private val aminoAcids: List<Char>) :
        NucleotideFragment(id, nucleotideLoci, nucleotideQuality, nucleotides, genes) {

    fun containsAminoAcid(index: Int): Boolean {
        return aminoAcidLoci.contains(index)
    }

    fun aminoAcid(loci: Int): Char {
        return aminoAcids[aminoAcidLoci.indexOf(loci)]
    }

    fun aminoAcidIndices(): List<Int> = aminoAcidLoci
}