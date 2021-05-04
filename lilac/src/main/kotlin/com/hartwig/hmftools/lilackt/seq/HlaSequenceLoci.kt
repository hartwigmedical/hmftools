package com.hartwig.hmftools.lilackt.seq

import com.hartwig.hmftools.lilackt.evidence.PhasedEvidence
import com.hartwig.hmftools.lilackt.hla.HlaAllele

data class HlaSequenceLoci(val allele: HlaAllele, val sequences: List<String>) {
    val length = sequences.size

    fun containsInserts(): Boolean {
        return sequences.any { it.length > 1 }
    }

    fun containsDeletes(): Boolean {
        return sequences.any { it == "." }
    }

    fun containsIndels(): Boolean {
        return sequences.any { it == "." || it.length > 1 }
    }

    fun sequence(locus: Int): String {
        return sequences[locus]
    }

    fun sequence(): String {
        return sequences.joinToString(separator = "").replace(".", "")
    }

    fun sequence(startLocus: Int, endLocus: Int): String {
        return sequences
                .filterIndexed { index, _ -> index in startLocus..endLocus }
                .joinToString("")
    }

    private fun sequence(vararg indices: Int): String {
        return indices.joinToString("") { if (it < sequences.size) sequences[it] else "*" }
    }

    override fun toString(): String {
        return "HlaSequenceLoci(allele=$allele, sequence=${sequence()})"
    }

    fun consistentWith(evidence: PhasedEvidence): Boolean {
        return consistentWithAny(evidence.evidence.keys, *evidence.aminoAcidIndices)
    }

    fun consistentWithAny(targetSequence: Collection<String>, vararg targetIndices: Int): Boolean {
        return targetSequence.any { consistentWith(it, *targetIndices) }
    }

    fun consistentWith(targetSequence: String, vararg targetIndices: Int): Boolean {
        return match(targetSequence, *targetIndices) != HlaSequenceMatch.NONE
    }

    fun match(targetSequence: String, vararg targetIndices: Int): HlaSequenceMatch {
        if (targetIndices.isEmpty()) {
            return HlaSequenceMatch.NONE
        }

        val hlaSequence = sequence(*targetIndices)
        if (hlaSequence.length != targetSequence.length) {
            return HlaSequenceMatch.NONE
        }

        var wildCardCount = 0
        for (index in targetSequence.indices) {
            val target = targetSequence[index]
            if (hlaSequence[index] != '*' && hlaSequence[index] != target) {
                return HlaSequenceMatch.NONE
            }
            if (hlaSequence[index] == '*') {
                wildCardCount++
            }
        }

        if (wildCardCount > 0) {
            return if (wildCardCount == targetIndices.size) return HlaSequenceMatch.WILD else HlaSequenceMatch.PARTIAL
        }

        return HlaSequenceMatch.FULL

    }

    companion object {

        fun create(sequences: List<HlaSequence>): List<HlaSequenceLoci> {
            val reference = sequences[0].rawSequence
            return sequences.map { create(it.allele, it.rawSequence, reference) }
        }

        fun create(allele: HlaAllele, sequence: String, reference: String): HlaSequenceLoci {
            val sequences = mutableListOf<String>()

            fun isBaseIgnored(i: Int) = (sequence[i] == '.' && reference[i] == '.') || (sequence[i] == '|')
            fun isBaseInserted(i: Int) = sequence[i] != '.' && (i >= reference.length || reference[i] == '.')
            var insLength = 0

            for (i in sequence.indices) {
                val isInsert = isBaseInserted(i)
                val isIgnored = isBaseIgnored(i)

                if (insLength > 0 && !isInsert) {
                    val insert = sequence.substring(i - insLength, i)
                    sequences[sequences.size - 1] = sequences[sequences.size - 1] + insert
                    insLength = 0
                }

                if (isInsert) {
                    insLength++
                } else if (!isIgnored) {
                    val locusSequence = sequence[i].toString()
                    if (locusSequence == "-") {
                        sequences.add(reference[i].toString())
                    } else {
                        sequences.add(locusSequence)
                    }
                }
            }
            return HlaSequenceLoci(allele, sequences)
        }
    }

}