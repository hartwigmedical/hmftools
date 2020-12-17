package com.hartwig.hmftools.lilac.coverage

class ExonCoverage(val exonSequence: String) {
    val coverage = Array(exonSequence.length) { 0 }

    fun addCoverageFromStart(n: Int) {
        for (i in 0 until n) {
            coverage[i] += 1
        }
    }

    fun addCoverageFromEnd(n: Int) {
        for (i in coverage.size - n until coverage.size) {
            coverage[i] += 1
        }
    }

    fun addCompleteCoverage() {
        for (i in exonSequence.indices) {
            coverage[i] += 1
        }
    }

    fun isCovered(): Boolean {
        return minCoverage() > 0
    }

    fun minCoverage(): Int {
        return coverage.min()!!
    }

    fun maxCoverage(): Int {
        return coverage.max()!!
    }
}