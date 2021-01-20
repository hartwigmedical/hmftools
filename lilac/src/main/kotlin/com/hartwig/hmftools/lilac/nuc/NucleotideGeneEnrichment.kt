package com.hartwig.hmftools.lilac.nuc

class NucleotideGeneEnrichment(aBoundaries: Set<Int>, bBoundaries: Set<Int>, cBoundaries: Set<Int>) {

    private val aFilterB = ((aBoundaries subtract bBoundaries) union (bBoundaries subtract aBoundaries)).aminoAcidToNucleotides()
    private val aFilterC = ((aBoundaries subtract cBoundaries) union (cBoundaries subtract aBoundaries)).aminoAcidToNucleotides()
    private val bFilterA = ((bBoundaries subtract aBoundaries) union (aBoundaries subtract bBoundaries)).aminoAcidToNucleotides()
    private val bFilterC = ((bBoundaries subtract cBoundaries) union (cBoundaries subtract bBoundaries)).aminoAcidToNucleotides()
    private val cFilterA = ((cBoundaries subtract aBoundaries) union (aBoundaries subtract cBoundaries)).aminoAcidToNucleotides()
    private val cFilterB = ((cBoundaries subtract bBoundaries) union (bBoundaries subtract cBoundaries)).aminoAcidToNucleotides()

    fun enrich(fragments: List<NucleotideFragment>): List<NucleotideFragment> {
        return fragments.map { it.enrichGenes() }
    }

    private fun NucleotideFragment.enrichGenes(): NucleotideFragment {
        val genes = mutableSetOf<String>()
        if (this.matchToA()) {
            genes.add("HLA-A")
        }

        if (this.matchToB()) {
            genes.add("HLA-B")
        }

        if (this.matchToC()) {
            genes.add("HLA-C")
        }

        return NucleotideFragment(this.id, this.nucleotideIndices(), this.nucleotideQuality(), this.nucleotides(), genes)
    }

    private fun NucleotideFragment.matchToA(): Boolean {
        return this.genes.contains("HLA-A")
                || (this.genes.contains("HLA-B") && !aFilterB.any { this.containsNucleotide(it) })
                || (this.genes.contains("HLA-C") && !aFilterC.any { this.containsNucleotide(it) })
    }

    private fun NucleotideFragment.matchToB(): Boolean {
        return this.genes.contains("HLA-B")
                || (this.genes.contains("HLA-A") && !bFilterA.any { this.containsNucleotide(it) })
                || (this.genes.contains("HLA-C") && !bFilterC.any { this.containsNucleotide(it) })
    }

    private fun NucleotideFragment.matchToC(): Boolean {
        return this.genes.contains("HLA-C")
                || (this.genes.contains("HLA-A") && !cFilterA.any { this.containsNucleotide(it) })
                || (this.genes.contains("HLA-B") && !cFilterB.any { this.containsNucleotide(it) })
    }

    private fun Set<Int>.aminoAcidToNucleotides(): Set<Int> {
        return this.flatMap { listOf(3 * it, 3 * it + 1, 3 * it + 2) }.toSet()
    }


}