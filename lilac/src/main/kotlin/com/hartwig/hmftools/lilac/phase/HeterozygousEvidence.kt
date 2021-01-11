package com.hartwig.hmftools.lilac.phase

import com.hartwig.hmftools.lilac.read.Fragment
import java.util.*

class HeterozygousEvidence(val minBaseQual: Int, val heterozygousIndices: List<Int>, val readFragments: List<Fragment>) {


    fun initialEvidence(): List<PhasedEvidence> {
        val result = mutableListOf<PhasedEvidence>()

        for (i in 0..(heterozygousIndices.size - 2) ) {

            val evidence = PhasedEvidence.evidence(minBaseQual, readFragments,
                    heterozygousIndices[i],
                    heterozygousIndices[i + 1]
//                    heterozygousIndices[i + 2]
//                    heterozygousIndices[i + 3]
//                    heterozygousIndices[i + 4],
//                    heterozygousIndices[i + 5]
            )
            if (evidence.evidence.isNotEmpty()) {
                result.add(evidence)
            }
        }

        return result.sorted().filter { it.totalEvidence() > 10 }
    }


    fun extendEvidence3(current: PhasedEvidence, others: Set<PhasedEvidence>): List<PhasedEvidence> {
        val existingIndices = current.aminoAcidIndices
        val remainingIndices = readFragments
                .flatMap { it.aminoAcidIndices() }
                .intersect(heterozygousIndices)
                .filter { it !in existingIndices }

        val result = mutableListOf<PhasedEvidence>()
        for (i in remainingIndices) {
            val newIndices = (existingIndices + i).sortedArray()
            val fake = PhasedEvidence(newIndices, Collections.emptyMap())
            if (!others.contains(fake)) {
                val evidence = PhasedEvidence.evidence(minBaseQual, readFragments, *newIndices)
                if (evidence.evidence.isNotEmpty() && evidence.minEvidence() > 2) {
                    result.add(evidence)
                }
            }
        }
        return result.sorted().take(1).filter { it.totalEvidence() > 20 }
    }

    fun extendEvidence(existingEvidence: PhasedEvidence, others: Set<PhasedEvidence>): List<PhasedEvidence> {
        val existingIndices = existingEvidence.aminoAcidIndices
        val remainingIndices = readFragments
                .flatMap { it.aminoAcidIndices() }
                .intersect(heterozygousIndices)
                .filter { it !in existingIndices }

        val minExisting = existingIndices.min()!!
        val maxExisting = existingIndices.max()!!

        val remainingIndicesAbove = remainingIndices.filter { it > maxExisting }.sorted()
        val remainingIndicesBelow = remainingIndices.filter { it < minExisting }.sorted().reversed()

        val result = mutableListOf<PhasedEvidence>()
        if (remainingIndicesAbove.isNotEmpty()) {
            val newIndices = existingIndices + remainingIndicesAbove[0]
            val fake = PhasedEvidence(newIndices, Collections.emptyMap())
            if (!others.contains(fake)) {
                val evidence = PhasedEvidence.evidence(minBaseQual, readFragments, *newIndices)
                if (evidence.evidence.isNotEmpty() && evidence.minEvidence() >= 3) {
                    result.add(evidence)
                }
            }
        }

        if (remainingIndicesBelow.isNotEmpty()) {
            val newIndices = (existingIndices + remainingIndicesBelow[0]).sortedArray()
            val fake = PhasedEvidence(newIndices, Collections.emptyMap())
            if (!others.contains(fake)) {

                val evidence = PhasedEvidence.evidence(minBaseQual, readFragments, *newIndices)
                if (evidence.evidence.isNotEmpty() && evidence.minEvidence() >= 3) {
                    result.add(evidence)
                }
            }
        }

        return result.sorted().filter { it.totalEvidence() > 20 }
    }

}