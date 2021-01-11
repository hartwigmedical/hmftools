package com.hartwig.hmftools.lilac.seq

import com.hartwig.hmftools.lilac.hla.HlaAllele
import com.hartwig.hmftools.lilac.phase.PhasedEvidence

data class HlaSequence(val contig: String, val rawSequence: String) {
    val allele = HlaAllele(contig)
    val sequence = rawSequence.aligned()
    val length = sequence.length

    fun consistentWith(evidence: PhasedEvidence): Boolean {
        return evidence.evidence.keys.any { this.consistentWith(evidence.aminoAcidIndices, it.toCharArray()) }
    }

    fun consistentWith(aminoAcidIndices: IntArray, aminoAcids: CharArray): Boolean {
        for (i in aminoAcidIndices.indices) {
            val index = aminoAcidIndices[i]
            val aminoAcid = aminoAcids[i]

            if (this.length > index && this.sequence[index] != '*' && this.sequence[index] != aminoAcid) {
                return false
            }
        }

        return true
    }

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
        for (i in sequence.indices) {
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
        return "${contig.padEnd(20, ' ')}\t${sequence}"
    }

    private fun String.aligned(): String {
        return this.replace(".", "").replace("|", "")
    }

    private fun String.removeSpecialCharacters(): String {
        return this.replace("*", "").replace(".", "").replace("|", "")
    }

}