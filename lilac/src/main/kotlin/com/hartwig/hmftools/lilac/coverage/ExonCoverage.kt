package com.hartwig.hmftools.lilac.coverage

class ExonCoverage(val exonSequence: String) {
    val coverage = Array(exonSequence.length) { 0 }

    fun length(): Int {
        return exonSequence.length
    }

    fun addCoverageFromStart(n: Int) {
        synchronized(coverage) {
            for (i in 0 until n) {
                coverage[i] += 1
            }
        }
    }

    fun addCoverageFromEnd(n: Int) {
        synchronized(coverage) {
            for (i in coverage.size - n until coverage.size) {
                coverage[i] += 1
            }
        }
    }

    fun addCompleteCoverage() {
        synchronized(coverage) {
            for (i in coverage.indices) {
                coverage[i] += 1
            }
        }
    }

    fun addCoverage(index: Int, length: Int) {
        synchronized(coverage) {
            for (i in 0 until length) {
                coverage[index + i] += 1
            }
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