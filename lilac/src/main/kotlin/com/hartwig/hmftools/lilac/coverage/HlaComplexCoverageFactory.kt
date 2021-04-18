package com.hartwig.hmftools.lilac.coverage

import com.hartwig.hmftools.common.progress.FutureProgressTracker
import com.hartwig.hmftools.lilac.LilacConfig
import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.read.FragmentAlleles
import com.hartwig.hmftools.lilac.read.FragmentAlleles.Companion.filter
import java.util.concurrent.Callable
import java.util.concurrent.ExecutorService
import java.util.concurrent.Future

class HlaComplexCoverageFactory(private val config: LilacConfig) {
    private val progressTracker = FutureProgressTracker(0.1, 10000)

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
        val topRankedKeepers = topRanked.filter { it in config.commonAlleles }

        return (topTakers + topRankedKeepers).distinct()
    }

    fun rankedComplexCoverage(executorService: ExecutorService, fragmentAlleles: List<FragmentAlleles>, complexes: List<HlaComplex>,  recovered: List<HlaAllele>): List<HlaComplexCoverage> {
        val ranking = HlaComplexCoverageRanking(config.maxDistanceFromTopScore, config.commonAlleles, recovered)

        val list = mutableListOf<Future<HlaComplexCoverage>>()
        for (complex in complexes) {
            val untrackedCallable: Callable<HlaComplexCoverage> = Callable { proteinCoverage(fragmentAlleles, complex.alleles) }
            val trackedCallable: Callable<HlaComplexCoverage> = progressTracker.add(untrackedCallable)
            list.add(executorService.submit(trackedCallable))
        }

        val result = list.map { it.get() }
        return ranking.candidateRanking(result)

    }


}