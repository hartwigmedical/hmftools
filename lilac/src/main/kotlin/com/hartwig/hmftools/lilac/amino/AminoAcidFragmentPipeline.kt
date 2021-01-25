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
        val spliceEnricher = NucleotideSpliceEnrichment(minBaseQuality, minBaseCount, aBoundaries)

        val geneSpecific = geneEnriched.filter { it.genes.contains("HLA-A") }
        val spliceEnriched = spliceEnricher.enrich(geneSpecific)
        val result =  heterozygousIntersector.heterozygousIntersect(spliceEnriched)

        val jon = result.filter { it.containsAminoAcid(3) && it.aminoAcid(3) == 'R'}

        return result
    }

    fun typeB(): List<AminoAcidFragment> {
        val spliceEnricher = NucleotideSpliceEnrichment(minBaseQuality, minBaseCount, bBoundaries)

        val geneSpecific = geneEnriched.filter { it.genes.contains("HLA-B") }
        val spliceEnriched = spliceEnricher.enrich(geneSpecific)
        return heterozygousIntersector.heterozygousIntersect(spliceEnriched)
    }

    fun typeC(): List<AminoAcidFragment> {
        val spliceEnricher = NucleotideSpliceEnrichment(minBaseQuality, minBaseCount, cBoundaries)

        val geneSpecific = geneEnriched.filter { it.genes.contains("HLA-C") }
        val spliceEnriched = spliceEnricher.enrich(geneSpecific)
        return heterozygousIntersector.heterozygousIntersect(spliceEnriched)
    }

    fun combined(): List<AminoAcidFragment> {
        val combinedBoundaries = (aBoundaries + bBoundaries + cBoundaries)
        val spliceEnricher = NucleotideSpliceEnrichment(minBaseQuality, minBaseCount, combinedBoundaries)
        val spliceEnriched = spliceEnricher.enrich(geneEnriched)
        return heterozygousIntersector.heterozygousIntersect(spliceEnriched)
    }

}