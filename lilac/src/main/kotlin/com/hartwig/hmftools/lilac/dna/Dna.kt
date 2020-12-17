package com.hartwig.hmftools.lilac.dna

import com.hartwig.hmftools.common.codon.Codons


fun String.dnaReverseComplement(): String {
    val builder = StringBuilder()
    for (i in this.indices.reversed()) {
        when (this[i]) {
            'G' -> builder.append('C')
            'C' -> builder.append('G')
            'A' -> builder.append('T')
            'T' -> builder.append('A')
            else -> builder.append(this[i])
        }
    }
    return builder.toString()
}

fun String.aminoAcids(): String {
    return Codons.aminoAcids(this)
}