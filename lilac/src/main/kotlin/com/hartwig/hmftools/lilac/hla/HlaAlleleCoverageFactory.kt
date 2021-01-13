package com.hartwig.hmftools.lilac.hla

import com.hartwig.hmftools.lilac.read.Fragment
import com.hartwig.hmftools.lilac.read.FragmentAlleles
import com.hartwig.hmftools.lilac.seq.HlaSequence

class HlaAlleleCoverageFactory(private val hetLoci: Collection<Int>, private val fragments: List<Fragment>) {

    companion object {
        fun List<HlaAlleleCoverage>.combinedCoverage(): Double {
            return this.map { it.combinedCoverage }.sum()
        }

        fun List<HlaAlleleCoverage>.uniqueCoverage(): Int {
            return this.map { it.uniqueCoverage }.sum()
        }

        fun List<HlaAlleleCoverage>.totalCoverage(): Double {
            return this.map { it.combinedCoverage + it.uniqueCoverage }.sum()
        }

    }


    fun alleleCoverage(alleles: Collection<HlaSequence>): List<HlaAlleleCoverage> {
        val fragmentAlleles = FragmentAlleles.create(fragments, hetLoci, alleles)
        return HlaAlleleCoverage.proteinCoverage(fragmentAlleles)
    }


}