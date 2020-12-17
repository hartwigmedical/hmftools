package com.hartwig.hmftools.lilac.coverage

import com.hartwig.hmftools.lilac.hla.HlaAllele

class ProteinCoverage(val hlaAllele: HlaAllele, proteins: List<String>) {
    val map = proteins.map { Pair(it, ExonCoverage(it)) }.toMap()

    fun isCovered(): Boolean {
        return map.all { it.value.isCovered() }
    }

    fun minCoverage(): Int {
        return map.map { it.value.minCoverage() }.min()!!
    }

    fun maxCoverage(): Int {
        return map.map { it.value.maxCoverage() }.max()!!
    }

}