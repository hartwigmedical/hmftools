package com.hartwig.hmftools.lilac.coverage

import com.hartwig.hmftools.common.progress.FutureProgressTracker
import com.hartwig.hmftools.lilac.amino.AminoAcidFragment
import com.hartwig.hmftools.lilac.hla.HlaAllele
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
    private val fragmentAlleles = FragmentAlleles.create(aminoAcidFragments, aminoAcidLoci, aminoAcidSequences, nucleotideLoci, nucleotideSequences)
    private val progressTracker = FutureProgressTracker(0.1, 10000)

    fun complexCoverage(complexes: List<HlaComplex>): List<HlaComplexCoverage> {
        val list = mutableListOf<Future<HlaComplexCoverage>>()
        for (complex in complexes) {
            val untrackedCallable: Callable<HlaComplexCoverage> = Callable { proteinCoverage(complex.alleles) }
            val trackedCallable: Callable<HlaComplexCoverage> = progressTracker.add(untrackedCallable)
            list.add(executorService.submit(trackedCallable))
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
        return FragmentAlleles.create(fragmentAlleles, alleles)
    }

}