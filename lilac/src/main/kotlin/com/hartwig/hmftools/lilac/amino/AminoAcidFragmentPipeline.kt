package com.hartwig.hmftools.lilac.amino

import com.hartwig.hmftools.lilac.LilacConfig
import com.hartwig.hmftools.lilac.hla.HlaContext
import com.hartwig.hmftools.lilac.nuc.NucleotideFragment
import com.hartwig.hmftools.lilac.nuc.NucleotideQualEnrichment
import com.hartwig.hmftools.lilac.nuc.NucleotideSpliceEnrichment

class AminoAcidFragmentPipeline(config: LilacConfig, private val geneEnriched: List<NucleotideFragment>) {
    private val minBaseQuality= config.minBaseQual
    private val minBaseCount = config.minEvidence
    private val aminoAcidEnricher = AminoAcidQualEnrichment(minBaseCount)
    private val nucleotideQualEnrichment = NucleotideQualEnrichment(minBaseQuality, minBaseCount)

    fun phasingFragments(context: HlaContext): List<AminoAcidFragment> {
        val gene = "HLA-${context.gene}"
        val geneSpecific = geneEnriched.filter { it.genes.contains(gene) }

        return process(context.aminoAcidBoundaries, geneSpecific)
    }

    fun coverageFragments(): List<AminoAcidFragment> {
        return qualFiltered(geneEnriched)
    }

    private fun process(boundaries: Set<Int>, fragments: List<NucleotideFragment>): List<AminoAcidFragment> {
        val spliceEnricher = NucleotideSpliceEnrichment(minBaseQuality, minBaseCount, boundaries)

        val qualEnriched = nucleotideQualEnrichment.enrich(fragments)
        val spliceEnriched = spliceEnricher.enrich(qualEnriched)
        val result = aminoAcidEnricher.enrich(spliceEnriched)

        return result
    }

    private fun qualFiltered(fragments: List<NucleotideFragment>): List<AminoAcidFragment> {
        val qualityFilteredFragments = fragments.map { it.qualityFilter(minBaseQuality) }.filter { it.isNotEmpty() }
        return qualityFilteredFragments.map { it.toAminoAcidFragment() }
    }

}