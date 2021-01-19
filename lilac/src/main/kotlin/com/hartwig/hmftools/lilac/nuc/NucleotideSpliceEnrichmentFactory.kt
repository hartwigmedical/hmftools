package com.hartwig.hmftools.lilac.nuc

class NucleotideSpliceEnrichmentFactory(private val minBaseCount: Int, private val aBoundaries: Set<Int>, private val bBoundaries: Set<Int>, private val cBoundaries: Set<Int>) {

    fun allNucleotides(nucleotideFragments: List<NucleotideFragment>): List<NucleotideFragment> {
        val allProteinExonBoundaries = (aBoundaries + bBoundaries + cBoundaries)
        return allEvidence(allProteinExonBoundaries, nucleotideFragments)
    }

    fun typeANucleotides(nucleotideFragments: List<NucleotideFragment>): List<NucleotideFragment> {
        return allEvidence(aBoundaries, nucleotideFragments.filter { it.genes.contains("HLA-A") })
    }

    fun typeBNucleotides(nucleotideFragments: List<NucleotideFragment>): List<NucleotideFragment> {
        return allEvidence(bBoundaries, nucleotideFragments.filter { it.genes.contains("HLA-B") })
    }

    fun typeCNucleotides(nucleotideFragments: List<NucleotideFragment>): List<NucleotideFragment> {
        return allEvidence(cBoundaries, nucleotideFragments.filter { it.genes.contains("HLA-C") })
    }


    private fun allEvidence(aminoAcidBoundaries: Set<Int>, nucleotideFragments: List<NucleotideFragment>): List<NucleotideFragment> {
        return NucleotideSpliceEnrichment(minBaseCount, aminoAcidBoundaries).enrich(nucleotideFragments)

    }
}