package com.hartwig.hmftools.lilac.hla

import com.hartwig.hmftools.lilac.read.Fragment
import com.hartwig.hmftools.lilac.read.FragmentAlleles
import com.hartwig.hmftools.lilac.seq.HlaSequence
import kotlin.math.roundToInt

class HlaAlleleCoverageFactory(
        private val fragments: List<Fragment>,
        private val aminoAcidLoci: Collection<Int>, private val aminoAcidSequences: Collection<HlaSequence>,
        private val nucleotideLoci: Collection<Int>, private val nucleotideSequences: Collection<HlaSequence>) {

    companion object {
        fun List<HlaAlleleCoverage>.combinedCoverage(): Double {
            return this.map { it.sharedCoverage }.sum()
        }

        fun List<HlaAlleleCoverage>.uniqueCoverage(): Int {
            return this.map { it.uniqueCoverage }.sum()
        }

        fun List<HlaAlleleCoverage>.totalCoverage(): Double {
            return this.map { it.sharedCoverage + it.uniqueCoverage }.sum()
        }

        fun List<HlaAlleleCoverage>.coverageString(confirmed: List<HlaAllele>): String {
            var shared = 0.0
            var unique = 0
            for (coverage in this) {
//                if (coverage.allele !in confirmed) {
                    shared += coverage.sharedCoverage
                    unique += coverage.uniqueCoverage
//                }
            }

            return "${(shared + unique).roundToInt()}\t${unique}\t${shared.roundToInt()}\t${this.sortedBy { x -> x.allele }.joinToString ("\t")}"

        }

    }

    fun groupCoverage(alleles: Collection<HlaAllele>): List<HlaAlleleCoverage> {
        val fragmentAlleles = fragmentAlleles(alleles)
        return HlaAlleleCoverage.groupCoverage(fragmentAlleles)
    }

    fun proteinCoverage(alleles: Collection<HlaAllele>): List<HlaAlleleCoverage> {
        val fragmentAlleles = fragmentAlleles(alleles)
        return HlaAlleleCoverage.proteinCoverage(fragmentAlleles)
    }

    private fun fragmentAlleles(alleles: Collection<HlaAllele>): List<FragmentAlleles> {
        val specificProteins = alleles.map { it.specificProtein() }
        val aminoAcids = aminoAcidSequences.filter { it.allele in alleles }
        val nucleotides = nucleotideSequences.filter { it.allele.specificProtein() in specificProteins }

        return FragmentAlleles.create(fragments, aminoAcidLoci, aminoAcids, nucleotideLoci, nucleotides)
    }


}