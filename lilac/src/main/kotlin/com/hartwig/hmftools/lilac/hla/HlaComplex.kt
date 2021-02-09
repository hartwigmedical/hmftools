package com.hartwig.hmftools.lilac.hla

data class HlaComplex(val alleles: List<HlaAllele>) {

    companion object {

        fun complexes(confirmedGroups: List<HlaAllele>, confirmedProteins: List<HlaAllele>, candidates: List<HlaAllele>): List<HlaComplex> {
            val a = gene(confirmedGroups.filter { it.gene == "A" }.take(2), confirmedProteins.filter { it.gene == "A" }.take(2), candidates.filter { it.gene == "A" })
            val b = gene(confirmedGroups.filter { it.gene == "B" }.take(2), confirmedProteins.filter { it.gene == "B" }.take(2), candidates.filter { it.gene == "B" })
            val c = gene(confirmedGroups.filter { it.gene == "C" }.take(2), confirmedProteins.filter { it.gene == "C" }.take(2), candidates.filter { it.gene == "C" })

            return combineComplexes(combineComplexes(a, b), c)
        }


        fun gene(confirmedGroups: List<HlaAllele>, confirmedProteins: List<HlaAllele>, candidates: List<HlaAllele>): List<HlaComplex> {

            if (confirmedProteins.size == 2) {
                return listOf(HlaComplex(confirmedProteins))
            }

            if (confirmedProteins.size == 1) {
                val confirmedProteinGroups = confirmedProteins.map { it.asAlleleGroup() }
                val remainingGroups = confirmedGroups.filter { it !in confirmedProteinGroups }
                val first = confirmedProteins
                val second = if (remainingGroups.isEmpty()) candidates.filter { it != confirmedProteins[0] } else candidates.filter { it.asAlleleGroup() in remainingGroups }

                return if (remainingGroups.isEmpty())  combineAlleles(first, second) + HlaComplex(first) else  combineAlleles(first, second)
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


