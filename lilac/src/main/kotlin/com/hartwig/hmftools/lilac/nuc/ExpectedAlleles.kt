package com.hartwig.hmftools.lilac.nuc

class ExpectedAlleles(private val counts: List<Pair<Int, Int>>) {

    fun expectedAlleles(loci: Collection<Int>): Int {
        return (counts.filter { it.first in loci }.map { it.second }.min() ?: 3) * 2
    }

    fun expectedAlleles(loci: Int): Int {
        return (counts.filter { it.first == loci }.map { it.second }.min() ?: 3) * 2
    }

    companion object {
        fun create(target: Set<Int>, other1: Set<Int>, other2: Set<Int>): ExpectedAlleles {
            val allLocations = (target union other1 union other2).sorted()
            return ExpectedAlleles(createCounts(allLocations, target, other1, other2))
        }

        private fun createCounts(allLocations: Collection<Int>, target: Set<Int>, other1: Set<Int>, other2: Set<Int>): List<Pair<Int, Int>> {
            val result = mutableListOf<Pair<Int, Int>>()
            for (loci in allLocations) {
                val inTarget = target.contains(loci)
                val inOther1 = other1.contains(loci).toInt()
                val inOther2 = other2.contains(loci).toInt()
                if (inTarget) {
                    result.add(Pair(loci, 1 + inOther1 + inOther2))
                } else {
                    result.add(Pair(loci, 3 - inOther1 - inOther2))
                }

            }

            return result
        }

        private fun Boolean.toInt(): Int {
            return if (this) 1 else 0
        }

    }
}


