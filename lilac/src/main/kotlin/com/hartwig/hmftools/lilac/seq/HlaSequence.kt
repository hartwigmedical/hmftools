package com.hartwig.hmftools.lilac.seq

import com.hartwig.hmftools.lilac.evidence.PhasedEvidence
import com.hartwig.hmftools.lilac.hla.HlaAllele

data class HlaSequence(val contig: String, val rawSequence: String) {
    val allele = HlaAllele(contig)
    val sequence = rawSequence.aligned()
    val length = sequence.length

    fun deletes(reference: String): List<HlaSequenceIndel> {
        val result = mutableListOf<HlaSequenceIndel>()

        var i = 0
        var delLoci = 0
        var delLength = 0
        fun inDelete(): Boolean = delLength != 0
        fun isBaseDeleted(loci: Int) = rawSequence[loci] == '.' && reference[loci] != '.'

        while (i < rawSequence.length) {
            val baseDeleted = isBaseDeleted(i)

            if (inDelete()) {
                if (baseDeleted) {
                    delLength++
                } else {
                    result.add(HlaSequenceIndel.HlaSequenceDelete(delLoci, delLength))
                    delLoci = 0
                    delLength = 0
                }
            } else if (baseDeleted){
                delLoci = i
                delLength = 1
            }

            i++
        }

        return result
    }

    fun inserts(reference: String): List<HlaSequenceIndel> {
        val result = mutableListOf<HlaSequenceIndel>()

        var i = 0
        var insertLoci = 0
        var insLength = 0
        fun inInsert(): Boolean = insLength != 0
        fun isBaseInserted(loci: Int) = rawSequence[loci] != '.' && (loci >= reference.length || reference[loci] == '.')

        while (i < rawSequence.length) {
            val baseInserted = isBaseInserted(i)

            if (inInsert()) {
                if (baseInserted) {
                    insLength++
                } else {
                    result.add(HlaSequenceIndel.HlaSequenceInsert(insertLoci, insLength))
                    insertLoci = 0
                    insLength = 0
                }
            } else if (baseInserted){
                insertLoci = i
                insLength = 1
            }

            i++
        }

        return result
    }

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