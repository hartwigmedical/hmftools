package com.hartwig.hmftools.lilac.amino

import com.hartwig.hmftools.lilac.hla.HlaContext
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment
import com.hartwig.hmftools.lilac.nuc.NucleotideSpliceEnrichment

class AminoAcidFragmentPipeline(private val minBaseQuality: Int, private val minBaseCount: Int, private val geneEnriched: List<NucleotideFragment>) {
    private val heterozygousIntersector = AminoAcidFragmentEnrichment(minBaseQuality, minBaseCount)

    fun type(context: HlaContext): List<AminoAcidFragment> {
        val gene = "HLA-${context.gene}"
        val spliceEnricher = NucleotideSpliceEnrichment(minBaseQuality, minBaseCount, context.aminoAcidBoundaries)

        val geneSpecific = geneEnriched.filter { it.genes.contains(gene) }
        val spliceEnriched = spliceEnricher.enrich(geneSpecific)
        val result = heterozygousIntersector.heterozygousIntersect(gene, spliceEnriched)

        return result
                .map { it.qualityFilterNucleotides(minBaseQuality) }
    }


    fun combined(combinedBoundaries: Set<Int>): List<AminoAcidFragment> {
        val spliceEnricher = NucleotideSpliceEnrichment(minBaseQuality, minBaseCount, combinedBoundaries)
        val spliceEnriched = spliceEnricher.enrich(geneEnriched)
        return heterozygousIntersector.heterozygousIntersect("ALL", spliceEnriched)
                .map { it.qualityFilterNucleotides(minBaseQuality) }
    }
}