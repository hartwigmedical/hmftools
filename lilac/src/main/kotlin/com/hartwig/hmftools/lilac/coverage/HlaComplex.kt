package com.hartwig.hmftools.lilac.coverage

import com.hartwig.hmftools.lilac.hla.HlaAllele

data class HlaComplex(val alleles: List<HlaAllele>) {

    companion object {

        fun complexes(confirmedGroups: List<HlaAllele>, confirmedProteins: List<HlaAllele>, candidates: List<HlaAllele>): List<HlaComplex> {
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

    }
}


