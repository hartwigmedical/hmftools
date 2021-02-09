package com.hartwig.hmftools.lilac.hla

import com.hartwig.hmftools.lilac.amino.AminoAcidFragment
import com.hartwig.hmftools.lilac.read.FragmentAlleles
import com.hartwig.hmftools.lilac.seq.HlaSequenceLoci
import java.util.concurrent.Callable
import java.util.concurrent.ExecutorService
import java.util.concurrent.Future

class HlaComplexCoverageFactory(
        private val executorService: ExecutorService,
        private val aminoAcidFragments: List<AminoAcidFragment>,
        private val aminoAcidLoci: Collection<Int>, private val aminoAcidSequences: Collection<HlaSequenceLoci>,
        private val nucleotideLoci: Collection<Int>, private val nucleotideSequences: Collection<HlaSequenceLoci>) {

    fun complexCoverage(complexes: List<HlaComplex>): List<HlaComplexCoverage> {
        val list = mutableListOf<Future<HlaComplexCoverage>>()
        for (complex in complexes) {
            val callable: Callable<HlaComplexCoverage> = Callable {
                proteinCoverage(complex.alleles)
            }
            list.add(executorService.submit(callable))
        }

        return list.map { it.get() }.sortedDescending()
    }

    fun groupCoverage(alleles: Collection<HlaAllele>): HlaComplexCoverage {
        val fragmentAlleles = fragmentAlleles(alleles)
        return HlaComplexCoverage.create(HlaAlleleCoverage.groupCoverage(fragmentAlleles))
    }

    fun proteinCoverage(alleles: Collection<HlaAllele>): HlaComplexCoverage {
        val fragmentAlleles = fragmentAlleles(alleles)
        return HlaComplexCoverage.create(HlaAlleleCoverage.proteinCoverage(fragmentAlleles))
    }

    private fun fragmentAlleles(alleles: Collection<HlaAllele>): List<FragmentAlleles> {
        val specificProteins = alleles.map { it.asFourDigit() }
        val aminoAcids = aminoAcidSequences.filter { it.allele in alleles }
        val nucleotides = nucleotideSequences.filter { it.allele.asFourDigit() in specificProteins }

        return FragmentAlleles.create(aminoAcidFragments, aminoAcidLoci, aminoAcids, nucleotideLoci, nucleotides)
    }


}