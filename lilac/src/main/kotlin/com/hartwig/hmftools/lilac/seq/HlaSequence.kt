package com.hartwig.hmftools.lilac.seq

import com.hartwig.hmftools.lilac.hla.HlaAllele

data class HlaSequence(val contig: String, val rawSequence: String) {
    val allele: HlaAllele by lazy { HlaAllele(contig) }
    val sequence = rawSequence.removeSpecialCharacters()
    val length = sequence.length

    fun copyWithAdditionalSequence(additionalSequence: String): HlaSequence {
        return HlaSequence(contig, rawSequence + additionalSequence)
    }

    fun inflate(template: String): HlaSequence {
        val joiner = StringBuilder()
        for (i in rawSequence.indices) {
            if (rawSequence[i] == '-') {
                joiner.append(template[i])
            } else {
                joiner.append(rawSequence[i])
            }
        }
        return HlaSequence(contig, joiner.toString())
    }

    fun deflate(template: String): HlaSequence {
        val joiner = StringBuilder()
        for (i in rawSequence.indices) {
            when {
                rawSequence[i] == '.' -> joiner.append('.')
                rawSequence[i] == '|' -> joiner.append('|')
                rawSequence[i] == '*' -> joiner.append('*')
                i > template.length - 1 -> joiner.append(rawSequence[i])
                rawSequence[i] == template[i] -> joiner.append('-')
                else -> joiner.append(rawSequence[i])
            }
        }
        return HlaSequence(contig, joiner.toString())
    }

    override fun toString(): String {
        return "${contig.padEnd(20, ' ')}\t${rawSequence}"
    }

    private fun String.removeSpecialCharacters(): String {
        return this.replace("*", "").replace(".", "").replace("|", "")
    }
}