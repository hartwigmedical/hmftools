package com.hartwig.hmftools.lilackt.nuc

class NucleotideGeneEnrichment(aBoundaries: Set<Int>, bBoundaries: Set<Int>, cBoundaries: Set<Int>) {

    val aFilterB = ((aBoundaries subtract bBoundaries) union (bBoundaries subtract aBoundaries)).min()!!
    val aFilterC = ((aBoundaries subtract cBoundaries) union (cBoundaries subtract aBoundaries)).min()!!
    val bFilterA = ((bBoundaries subtract aBoundaries) union (aBoundaries subtract bBoundaries)).min()!!
    val bFilterC = ((bBoundaries subtract cBoundaries) union (cBoundaries subtract bBoundaries)).min()!!
    val cFilterA = ((cBoundaries subtract aBoundaries) union (aBoundaries subtract cBoundaries)).min()!!
    val cFilterB = ((cBoundaries subtract bBoundaries) union (bBoundaries subtract cBoundaries)).min()!!

    fun enrich(fragments: List<NucleotideFragment>): List<NucleotideFragment> {
        return fragments.map { it.enrichGenes() }
    }

    fun enrich(fragment: NucleotideFragment): NucleotideFragment {
        return fragment.enrichGenes()
    }

    fun NucleotideFragment.enrichGenes(): NucleotideFragment {
        if (this.containsIndel()) {
            return this
        }

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

        return NucleotideFragment(this.id, genes, this.nucleotideLoci(), this.nucleotideQuality(), this.nucleotides())
    }

    private fun NucleotideFragment.matchToA(): Boolean {
        return this.genes.contains("HLA-A")
                || (this.genes.contains("HLA-B") && (this.nucleotideLoci().max() ?: 0) < 3 * aFilterB)
                || (this.genes.contains("HLA-C") && (this.nucleotideLoci().max() ?: 0) < 3 * aFilterC)
    }

    private fun NucleotideFragment.matchToB(): Boolean {
        return this.genes.contains("HLA-B")
                || (this.genes.contains("HLA-A") && (this.nucleotideLoci().max() ?: 0) < 3 * bFilterA)
                || (this.genes.contains("HLA-C") && (this.nucleotideLoci().max() ?: 0) < 3 * bFilterC)
    }

    private fun NucleotideFragment.matchToC(): Boolean {
        return this.genes.contains("HLA-C")
                || (this.genes.contains("HLA-A") && (this.nucleotideLoci().max() ?: 0) < 3 * cFilterA)
                || (this.genes.contains("HLA-B") && (this.nucleotideLoci().max() ?: 0) < 3 * cFilterB)
    }
}