package com.hartwig.hmftools.lilackt.nuc

import junit.framework.Assert.assertTrue
import org.junit.Test

class NucleotideGeneEnrichmentTest {

    private val aProteinExonBoundaries = setOf(24, 114, 206, 298, 337, 348, 364, 365)
    private val bProteinExonBoundaries = setOf(24, 114, 206, 298, 337, 348, 362)
    private val cProteinExonBoundaries = setOf(24, 114, 206, 298, 338, 349, 365, 366)
    private val enrichment = NucleotideGeneEnrichment(aProteinExonBoundaries, bProteinExonBoundaries, cProteinExonBoundaries)

    @Test
    fun testGeneEnrichment() {
        assertGene(setOf("HLA-A", "HLA-B"), "HLA-A", listOf(337, 348))
        assertGene(setOf("HLA-A", "HLA-B"), "HLA-B", listOf(337, 348))
        assertGene(setOf("HLA-C"), "HLA-C", listOf(337, 348))

        assertGene(setOf("HLA-A"), "HLA-A", listOf(337, 348, 362))
        assertGene(setOf("HLA-B"), "HLA-B", listOf(337, 348, 362))
    }

    fun assertGene(expectedGenes: Set<String>, alignedGene: String, aminoAcideIndices: List<Int>) {
        val victim = create(alignedGene, aminoAcideIndices.flatMap { listOf(it * 3, it * 3 + 1, it * 3 + 2) })
        val result = enrichment.enrich(listOf(victim))[0]
        assertTrue(result.genes == expectedGenes)
    }

    fun create(gene: String, indices: List<Int>): NucleotideFragment {
        return NucleotideFragment("id", setOf(gene), indices, indices.map { 0 }, indices.map { "G" })
    }


}