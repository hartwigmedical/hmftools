package com.hartwig.hmftools.lilackt.variant

import com.hartwig.hmftools.common.variant.VariantContextDecorator
import com.hartwig.hmftools.lilackt.LilacConfig
import com.hartwig.hmftools.lilackt.LociPosition
import com.hartwig.hmftools.lilac.coverage.HlaAlleleCoverage
import com.hartwig.hmftools.lilackt.read.FragmentAlleles
import com.hartwig.hmftools.lilackt.read.SAMRecordReader
import com.hartwig.hmftools.lilackt.seq.HlaSequenceLoci

class SomaticAlleleCoverage(private val config: LilacConfig, hetLoci: Collection<Int>, lociPosition: LociPosition, variants: List<VariantContextDecorator>, private val winners: Set<HlaSequenceLoci>) {

    private val variantLoci = variants
            .map { lociPosition.nucelotideLoci(it.position().toInt()) }
            .filter { it >= 0 }
            .map { it / 3 }

    private val hetLociSansVariants = hetLoci subtract variantLoci

    fun alleleCoverage(variant: VariantContextDecorator, reader: SAMRecordReader): List<HlaAlleleCoverage> {
        val variantFragments = reader
                .readFromBam(variant)
                .map { it.qualityFilter(config.minBaseQual) }
                .filter { it.isNotEmpty() }
                .map { it.toAminoAcidFragment() }


        val variantFragmentAlleles = FragmentAlleles.create(variantFragments, hetLociSansVariants, winners, listOf(), listOf())
        val coverage = HlaAlleleCoverage.proteinCoverage(variantFragmentAlleles).sortedByDescending { it.totalCoverage }
        return coverage.filter { it.totalCoverage == coverage[0].totalCoverage }
    }


}