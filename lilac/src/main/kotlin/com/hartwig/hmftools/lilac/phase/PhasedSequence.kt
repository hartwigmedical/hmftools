package com.hartwig.hmftools.lilac.phase

import com.hartwig.hmftools.lilac.read.SAMRecordRead

data class PhasedSequence(val startIndex: Int, val endIndex: Int, val evidence: Int, val sequence: String) {


    companion object {

        fun locallyPhased(startIndex: Int, endIndex: Int, alignments: List<SAMRecordRead>): List<PhasedSequence> {
            val locallyPhased = alignments
                    .filter { it.containsAminoAcid(startIndex) && it.containsAminoAcid(endIndex) }
                    .groupBy { it.aminoAcids(startIndex, endIndex, 30) }
                    .mapValues { (_, v) -> v.size }
                    .filter { it.key.isNotEmpty() }

            return locallyPhased.map { (sequence, count) -> PhasedSequence(startIndex, endIndex, count, sequence) }.sortedBy { -it.evidence }
        }
    }

}