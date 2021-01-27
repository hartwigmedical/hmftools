package com.hartwig.hmftools.lilac.amino

import com.hartwig.hmftools.lilac.nuc.NucleotideFragment
import com.hartwig.hmftools.lilac.nuc.NucleotideGeneEnrichment
import com.hartwig.hmftools.lilac.nuc.NucleotideSpliceEnrichment

class AminoAcidFragmentPipeline(private val minBaseQuality: Int, private val minBaseCount: Int,
                                private val aBoundaries: Set<Int>, private val bBoundaries: Set<Int>, private val cBoundaries: Set<Int>,
                                nucleotideFragment: List<NucleotideFragment>) {
    private val geneEnriched = NucleotideGeneEnrichment(aBoundaries, bBoundaries, cBoundaries).enrich(nucleotideFragment)
    private val heterozygousIntersector = AminoAcidFragmentEnrichment(minBaseQuality, minBaseCount)

    fun typeA(): List<AminoAcidFragment> {
        val gene = "HLA-A"
        val spliceEnricher = NucleotideSpliceEnrichment(minBaseQuality, minBaseCount, aBoundaries)

        val geneSpecific = geneEnriched.filter { it.genes.contains(gene) }
        val spliceEnriched = spliceEnricher.enrich(geneSpecific)
        val result =  heterozygousIntersector.heterozygousIntersect(gene, spliceEnriched)

        return result
                .map { it.qualityFilterNucleotides(minBaseQuality) }
    }

    fun typeB(): List<AminoAcidFragment> {
        val gene = "HLA-B"
        val spliceEnricher = NucleotideSpliceEnrichment(minBaseQuality, minBaseCount, bBoundaries)

        val geneSpecific = geneEnriched.filter { it.genes.contains(gene) }
        val spliceEnriched = spliceEnricher.enrich(geneSpecific)
        return heterozygousIntersector.heterozygousIntersect(gene, spliceEnriched)
                .map { it.qualityFilterNucleotides(minBaseQuality) }
    }

    fun typeC(): List<AminoAcidFragment> {
        val gene = "HLA-C"
        val spliceEnricher = NucleotideSpliceEnrichment(minBaseQuality, minBaseCount, cBoundaries)

        val geneSpecific = geneEnriched.filter { it.genes.contains(gene) }
        val spliceEnriched = spliceEnricher.enrich(geneSpecific)
        return heterozygousIntersector.heterozygousIntersect(gene, spliceEnriched)
                .map { it.qualityFilterNucleotides(minBaseQuality) }
    }

    fun combined(): List<AminoAcidFragment> {
        val combinedBoundaries = (aBoundaries + bBoundaries + cBoundaries)
        val spliceEnricher = NucleotideSpliceEnrichment(minBaseQuality, minBaseCount, combinedBoundaries)
        val spliceEnriched = spliceEnricher.enrich(geneEnriched)
        return heterozygousIntersector.heterozygousIntersect("ALL", spliceEnriched)
                .map { it.qualityFilterNucleotides(minBaseQuality) }
    }
}