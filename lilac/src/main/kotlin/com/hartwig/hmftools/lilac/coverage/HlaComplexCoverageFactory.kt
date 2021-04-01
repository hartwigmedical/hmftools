package com.hartwig.hmftools.lilac.coverage

import com.hartwig.hmftools.common.progress.FutureProgressTracker
import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.read.FragmentAlleles
import com.hartwig.hmftools.lilac.read.FragmentAlleles.Companion.filter
import java.util.concurrent.Callable
import java.util.concurrent.ExecutorService
import java.util.concurrent.Future

class HlaComplexCoverageFactory(maxDistanceFromTopScore: Int, private val common: List<HlaAllele>) {
    private val progressTracker = FutureProgressTracker(0.1, 10000)
    private val ranking = HlaComplexCoverageRanking(maxDistanceFromTopScore, common)

    companion object {
        fun groupCoverage(fragmentAlleles: List<FragmentAlleles>, alleles: Collection<HlaAllele>): HlaComplexCoverage {
            val filteredFragments = fragmentAlleles(fragmentAlleles, alleles)
            return HlaComplexCoverage.create(HlaAlleleCoverage.groupCoverage(filteredFragments))
        }

        fun proteinCoverage(fragmentAlleles: List<FragmentAlleles>, alleles: Collection<HlaAllele>): HlaComplexCoverage {
            val filteredFragments = fragmentAlleles(fragmentAlleles, alleles)
            return HlaComplexCoverage.create(HlaAlleleCoverage.proteinCoverage(filteredFragments))
        }

        private fun fragmentAlleles(fragmentAlleles: List<FragmentAlleles>, alleles: Collection<HlaAllele>): List<FragmentAlleles> {
            return fragmentAlleles.filter(alleles)
        }
    }

    fun rankedGroupCoverage(take: Int, fragmentAlleles: List<FragmentAlleles>, complexes: List<HlaComplex>): List<HlaAllele> {
        val topRanked =  complexes
                .map { proteinCoverage(fragmentAlleles, it.alleles) }
                .sortedBy { -it.totalCoverage }
                .flatMap { it.alleleCoverage }
                .map { it.allele }
                .distinct()

        val topTakers = topRanked.take(take)
        val topRankedKeepers = topRanked.filter { it in common }

        return (topTakers + topRankedKeepers).distinct()
    }

    fun rankedComplexCoverage(executorService: ExecutorService, fragmentAlleles: List<FragmentAlleles>, complexes: List<HlaComplex>): List<HlaComplexCoverage> {
        val list = mutableListOf<Future<HlaComplexCoverage>>()
        for (complex in complexes) {
            val untrackedCallable: Callable<HlaComplexCoverage> = Callable { proteinCoverage(fragmentAlleles, complex.alleles) }
            val trackedCallable: Callable<HlaComplexCoverage> = progressTracker.add(untrackedCallable)
            list.add(executorService.submit(trackedCallable))
        }

        val result = mutableListOf<HlaComplexCoverage>()
        val iterator = list.iterator()
        while (iterator.hasNext()) {
            val coverage = iterator.next().get()
            result.add(coverage)

            if (result.size > 100) {
                val filtered = ranking.candidateRanking(result)
                result.clear()
                result.addAll(filtered)
            }

            iterator.remove()
        }

        return ranking.candidateRanking(result)

    }


}