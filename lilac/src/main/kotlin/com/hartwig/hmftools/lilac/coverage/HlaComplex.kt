package com.hartwig.hmftools.lilac.coverage

import com.hartwig.hmftools.lilac.LilacConfig
import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.read.FragmentAlleles
import org.apache.logging.log4j.LogManager

data class HlaComplex(val alleles: List<HlaAllele>) {

    companion object {

        val logger = LogManager.getLogger(this::class.java)

        fun complexes(config: LilacConfig, referenceFragmentAlleles: List<FragmentAlleles>, candidateAlleles: List<HlaAllele>, recoveredAlleles: List<HlaAllele>): List<HlaComplex>  {

            logger.info("Identifying uniquely identifiable groups and proteins [total,unique,shared,wide]")
            val groupCoverage = HlaComplexCoverageFactory.groupCoverage(referenceFragmentAlleles, candidateAlleles)
            val confirmedGroups = groupCoverage.confirmUnique(config)
            val discardedGroups = groupCoverage.alleleCoverage.filter { it.uniqueCoverage > 0 && it !in confirmedGroups }.sortedDescending()
            if (confirmedGroups.isNotEmpty()) {
                logger.info("    confirmed ${confirmedGroups.size} unique groups: " + confirmedGroups.joinToString(", "))
            } else {
                logger.info("    confirmed 0 unique groups")
            }
            if (discardedGroups.isNotEmpty()) {
                logger.info("    found ${discardedGroups.size} insufficiently unique groups: " + discardedGroups.joinToString(", "))
            }

            val confirmedGroupAlleles = confirmedGroups.alleles()
            val candidatesAfterConfirmedGroups = candidateAlleles.filterWithConfirmedGroups(confirmedGroupAlleles)
            val proteinCoverage = HlaComplexCoverageFactory.proteinCoverage(referenceFragmentAlleles, candidatesAfterConfirmedGroups)
            val confirmedProtein = proteinCoverage.confirmUnique(config)
            val discardedProtein = proteinCoverage.alleleCoverage.filter { it.uniqueCoverage > 0 && it !in confirmedProtein }.sortedDescending()
            if (confirmedProtein.isNotEmpty()) {
                logger.info("    confirmed ${confirmedProtein.size} unique proteins: " + confirmedProtein.joinToString(", "))
            } else {
                logger.info("    confirmed 0 unique proteins")
            }
            if (discardedProtein.isNotEmpty()) {
                logger.info("    found ${discardedProtein.size} insufficiently unique proteins: " + discardedProtein.joinToString(", "))
            }

            val candidatesAfterConfirmedProteins = candidatesAfterConfirmedGroups.filterWithConfirmedProteins(confirmedProtein.map { it.allele })
            val aOnlyComplexes = gene("A", confirmedGroups.alleles(), confirmedProtein.alleles(), candidatesAfterConfirmedProteins)
            val bOnlyComplexes = gene("B", confirmedGroups.alleles(), confirmedProtein.alleles(), candidatesAfterConfirmedProteins)
            val cOnlyComplexes = gene("C", confirmedGroups.alleles(), confirmedProtein.alleles(), candidatesAfterConfirmedProteins)

            val complexes: List<HlaComplex>
            val simpleComplexCount = aOnlyComplexes.size.toLong() * bOnlyComplexes.size * cOnlyComplexes.size
            complexes = if (simpleComplexCount > 100_000 || simpleComplexCount < 0) {
                logger.info("Candidate permutations exceeds maximum complexity")
                val groupRankedCoverageFactory = HlaComplexCoverageFactory(config)
                val aTopCandidates = groupRankedCoverageFactory.rankedGroupCoverage(10, referenceFragmentAlleles, aOnlyComplexes)
                val bTopCandidates = groupRankedCoverageFactory.rankedGroupCoverage(10, referenceFragmentAlleles, bOnlyComplexes)
                val cTopCandidates = groupRankedCoverageFactory.rankedGroupCoverage(10, referenceFragmentAlleles, cOnlyComplexes)
                val topCandidates  = aTopCandidates + bTopCandidates + cTopCandidates
                val rejected = candidatesAfterConfirmedProteins subtract topCandidates
                logger.info("    discarding ${rejected.size} unlikely candidates: " + rejected.joinToString(", "))
                complexes(confirmedGroups.alleles(), confirmedProtein.alleles(), topCandidates)
            } else {
                complexes(confirmedGroups.alleles(), confirmedProtein.alleles(), candidatesAfterConfirmedProteins)
            }

            return complexes
        }

        private fun complexes(confirmedGroups: List<HlaAllele>, confirmedProteins: List<HlaAllele>, candidates: List<HlaAllele>): List<HlaComplex> {
            val a = gene("A", confirmedGroups, confirmedProteins, candidates)
            val b = gene("B", confirmedGroups, confirmedProteins, candidates)
            val c = gene("C", confirmedGroups, confirmedProteins, candidates)

            return combineComplexes(combineComplexes(a, b), c)
        }

        fun gene(gene: String, unfilteredGroups: List<HlaAllele>, unfilteredProteins: List<HlaAllele>, unfilteredCandidates: List<HlaAllele>): List<HlaComplex> {

            val confirmedGroups: List<HlaAllele> = unfilteredGroups.filter { it.gene == gene }.take(2)
            val confirmedProteins: List<HlaAllele> = unfilteredProteins.filter { it.gene == gene }.take(2)
            val candidates: List<HlaAllele> = unfilteredCandidates.filter { it.gene == gene }

            if (confirmedProteins.size == 2) {
                return listOf(HlaComplex(confirmedProteins))
            }

            if (confirmedProteins.size == 1) {
                val confirmedProteinGroups = confirmedProteins.map { it.asAlleleGroup() }
                val remainingGroups = confirmedGroups.filter { it !in confirmedProteinGroups }
                val first = confirmedProteins
                val second = if (remainingGroups.isEmpty()) candidates.filter { it != confirmedProteins[0] } else candidates.filter { it.asAlleleGroup() in remainingGroups }

                return if (remainingGroups.isEmpty()) combineAlleles(first, second) + HlaComplex(first) else combineAlleles(first, second)
            }

            if (confirmedGroups.size == 2) {
                val first = candidates.filter { it.asAlleleGroup() == confirmedGroups[0] }
                val second = candidates.filter { it.asAlleleGroup() == confirmedGroups[1] }
                return combineAlleles(first, second)
            }

            if (confirmedGroups.size == 1) {
                val first = candidates.filter { it.asAlleleGroup() == confirmedGroups[0] }
                val second = candidates
                return first.map { HlaComplex(listOf(it)) } + combineAlleles(first, second)
            }

            return candidates.map { HlaComplex(listOf(it)) } + combineAlleles(candidates, candidates)
        }

        fun combineComplexes(first: List<HlaComplex>, second: List<HlaComplex>): List<HlaComplex> {
            val intermediate = cartesianProduct(first, second)
            return intermediate.map { it.flatMap { it.alleles } }.map { HlaComplex(it) }
        }

        private fun combineAlleles(first: List<HlaAllele>, second: List<HlaAllele>): List<HlaComplex> {
            return cartesianProduct(first, second)
                    .map { it.sorted() }
                    .distinct()
                    .map { HlaComplex(it) }
        }

        private fun <T> cartesianProduct(first: Collection<T>, second: Collection<T>): List<List<T>> {
            val result = mutableListOf<List<T>>()
            for (i in first) {
                for (j in second) {
                    if (i != j) {
                        result.add(mutableListOf(i, j))
                    }
                }
            }
            return result.distinct()
        }

        private fun List<HlaAllele>.filterWithConfirmedProteins(confirmedGroups: List<HlaAllele>): List<HlaAllele> {
            return this.filterWithConfirmed(confirmedGroups) { it}
        }

        private fun List<HlaAllele>.filterWithConfirmedGroups(confirmedGroups: List<HlaAllele>): List<HlaAllele> {
            return this.filterWithConfirmed(confirmedGroups) { it.asAlleleGroup()}
        }

        private fun List<HlaAllele>.filterWithConfirmed(confirmed: List<HlaAllele>, transform: (HlaAllele) -> HlaAllele): List<HlaAllele> {
            val a = confirmed.filter { it.gene == "A" }
            val b = confirmed.filter { it.gene == "B" }
            val c = confirmed.filter { it.gene == "C" }
            val map = mapOf(Pair("A", a), Pair("B", b), Pair("C", c))

            return this.filter { map[it.gene]!!.size < 2 || map[it.gene]!!.contains(transform(it)) }
        }

        private fun HlaComplexCoverage.confirmUnique(config: LilacConfig): List<HlaAlleleCoverage> {
            val unique = this.alleleCoverage.filter { it.uniqueCoverage >= config.minConfirmedUniqueCoverage }.sortedDescending()
            val a = unique.filter { it.allele.gene == "A" }.take(2)
            val b = unique.filter { it.allele.gene == "B" }.take(2)
            val c = unique.filter { it.allele.gene == "C" }.take(2)

            return (a + b + c).sortedDescending()
        }

        private fun List<HlaAlleleCoverage>.alleles(): List<HlaAllele> {
            return this.map { it.allele }
        }
    }
}


