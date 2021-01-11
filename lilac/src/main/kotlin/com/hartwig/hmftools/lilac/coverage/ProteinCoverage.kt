package com.hartwig.hmftools.lilac.coverage

import com.hartwig.hmftools.lilac.hla.HlaAllele

class ProteinCoverage(val hlaAllele: HlaAllele, proteins: List<String>) {
    val map = proteins.map { Pair(it, ExonCoverage(it)) }.toMap() //TODO: REMOVE
    val coverages = proteins.map { ExonCoverage(it) }

    fun isCovered(): Boolean {
        return coverages.all { it.isCovered() }
    }

    fun minCoverage(): Int {
        return coverages.map { it.minCoverage() }.min()!!
    }

    fun maxCoverage(): Int {
        return coverages.map { it.maxCoverage() }.max()!!
    }

    fun avgCoverage(): Double {
        var total: Int = 0
        var sum: Int = 0

        for (coverage in coverages) {
            total += coverage.length()
            sum += coverage.coverage.sum()
        }

        return 1.0 * sum / total
    }

}