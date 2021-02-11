package com.hartwig.hmftools.lilac.amino

import com.hartwig.hmftools.lilac.hla.HlaContext
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment
import com.hartwig.hmftools.lilac.nuc.NucleotideQualEnrichment
import com.hartwig.hmftools.lilac.nuc.NucleotideSpliceEnrichment

class AminoAcidFragmentPipeline(private val minBaseQuality: Int, private val minBaseCount: Int, private val geneEnriched: List<NucleotideFragment>) {
    private val aminoAcidEnricher = AminoAcidQualEnrichment(minBaseCount)
    private val nucleotideQualEnrichment = NucleotideQualEnrichment(minBaseQuality, minBaseCount)

    fun phasingFragments(context: HlaContext): List<AminoAcidFragment> {
        val gene = "HLA-${context.gene}"
        val geneSpecific = geneEnriched.filter { it.genes.contains(gene) }

        return process(context.aminoAcidBoundaries, geneSpecific)
    }

    fun coverageFragments(combinedBoundaries: Set<Int>): List<AminoAcidFragment> {
        return process(combinedBoundaries, geneEnriched)
                .map { it.qualityFilterNucleotides(minBaseQuality) }
    }

    fun process(boundaries: Set<Int>, fragments: List<NucleotideFragment>): List<AminoAcidFragment> {
        val spliceEnricher = NucleotideSpliceEnrichment(minBaseQuality, minBaseCount, boundaries)

        val qualEnriched = nucleotideQualEnrichment.enrich(fragments)
        val spliceEnriched = spliceEnricher.enrich(qualEnriched)
        val result = aminoAcidEnricher.enrich(spliceEnriched)

        return result
    }

}