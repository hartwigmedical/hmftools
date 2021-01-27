package com.hartwig.hmftools.lilac.seq

import com.hartwig.hmftools.lilac.evidence.PhasedEvidence
import com.hartwig.hmftools.lilac.hla.HlaAllele

data class HlaSequence(val contig: String, val rawSequence: String) {
    val allele = HlaAllele(contig)
    val sequence = rawSequence.aligned()
    val length = sequence.length


    fun consistentWith(evidence: PhasedEvidence): Boolean {
        return evidence.evidence.keys.any { this.consistentWith(evidence.aminoAcidIndices, it.toCharArray()) }
    }

    fun consistentWith(sequenceIndices: IntArray, sequences: Collection<CharArray>): Boolean {
        return sequences.any { consistentWith(sequenceIndices, it) }
    }

    fun consistentWith(sequenceIndices: IntArray, sequence: CharArray): Boolean {
        return match(sequenceIndices, sequence) != HlaSequenceMatch.NONE
    }


    fun match(indicies: IntArray, sequence: CharArray): HlaSequenceMatch {
        if (indicies.isEmpty()) {
            return HlaSequenceMatch.NONE
        }

        var wildCardCount = 0
        for (i in indicies.indices) {
            val index = indicies[i]
            val aminoAcid = sequence[i]

            if (this.length > index && this.sequence[index] != '*' && this.sequence[index] != aminoAcid) {
                return HlaSequenceMatch.NONE
            }

            if (this.length > index && this.sequence[index] == '*') {
                wildCardCount++
            }
        }

        if (wildCardCount > 0) {
            return if (wildCardCount == indicies.size) return HlaSequenceMatch.WILD else HlaSequenceMatch.PARTIAL
        }

        return HlaSequenceMatch.FULL
    }

    fun copyWithAdditionalSequence(additionalSequence: String): HlaSequence {
        return HlaSequence(contig, rawSequence + additionalSequence)
    }

    fun pad(length: Int): HlaSequence {
        val currentLength = sequence.length
        val desiredLength = rawSequence.length + length - currentLength

        return HlaSequence(contig, rawSequence.padEnd(desiredLength, '*'))
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
                sequence[i] == '.' -> joiner.append('.')
                sequence[i] == '|' -> joiner.append('|')
                sequence[i] == '*' -> joiner.append('*')
                i > template.length - 1 -> joiner.append(sequence[i])
                sequence[i] == template[i] -> joiner.append('-')
                else -> joiner.append(sequence[i])
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