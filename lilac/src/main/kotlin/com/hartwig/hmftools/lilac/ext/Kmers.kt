package com.hartwig.hmftools.lilac.ext

import java.util.*

fun String.rollingKmers(minKmerLength: Int, maxKmerLength: Int): List<String> {
    if (this.length < minKmerLength) {
        return Collections.emptyList()
    }

    return this.rollingKmers(maxKmerLength.coerceAtMost(this.length))
}

fun String.rollingKmers(kmerSize: Int): List<String> {
    val result = mutableListOf<String>()
    for (i in 0..this.length - kmerSize) {
        result.add(this.substring(i, i + kmerSize))
    }

    return result
}
