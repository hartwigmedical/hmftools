package com.hartwig.hmftools.paddle.mutation

import com.hartwig.hmftools.paddle.Impact
import com.hartwig.hmftools.paddle.dnds.DndsMutation
import com.hartwig.hmftools.paddle.dnds.DndsMutationComparator
import org.apache.logging.log4j.LogManager

data class Mutations(val known: Int, val unknown: Int) {
    val total = known + unknown

    fun combine(other: Mutations): Mutations {
        return Mutations(known + other.known, unknown + other.unknown)
    }
}

data class MutationsGene(
        val gene: String, val synonymous: Int, val redundant: Int,
        val missense: Mutations, val nonsense: Mutations, val splice: Mutations, val inframe: Mutations, val frameshift: Mutations) {


    fun add(impact: Impact, known: Int, unknown: Int): MutationsGene {
        return when (impact) {
            Impact.MISSENSE -> copy(missense = missense.combine(Mutations(known, unknown)))
            Impact.NONSENSE -> copy(nonsense = nonsense.combine(Mutations(known, unknown)))
            Impact.SPLICE -> copy(splice = splice.combine(Mutations(known, unknown)))
            Impact.INFRAME -> copy(inframe = inframe.combine(Mutations(known, unknown)))
            Impact.FRAMESHIFT -> copy(frameshift = frameshift.combine(Mutations(known, unknown)))
            else -> throw IllegalStateException("Unexpected impact: $impact")
        }
    }

    fun combine(count: MutationsGene): MutationsGene {
        if (gene != count.gene) {
            throw IllegalArgumentException("Incompatible genes: $gene != ${count.gene}")
        }

        return MutationsGene(gene, synonymous + count.synonymous, redundant + count.redundant,
                missense.combine(count.missense), nonsense.combine(count.nonsense), splice.combine(count.splice), inframe.combine(count.inframe), frameshift.combine(count.frameshift)
        )
    }

    companion object {
        val logger = LogManager.getLogger(this::class.java)

        private val tsgComparator = DndsMutationComparator(true)
        private val oncoComparator = DndsMutationComparator(false)
        private val emptyCount = Mutations(0, 0)
        private val empty = MutationsGene("empty", 0, 0, emptyCount, emptyCount, emptyCount, emptyCount, emptyCount)

        fun oncoGeneMutations(mutations: List<DndsMutation>): List<MutationsGene> {
            return summary(mutations, this::oncoSampleSummary)
        }

        fun tsgGeneMutations(mutations: List<DndsMutation>): List<MutationsGene> {
            return summary(mutations, this::tsgSampleSummary)
        }

        private fun oncoSampleSummary(gene: String, sampleMutations: List<DndsMutation>): MutationsGene {
            fun Boolean.toInt() = if (this) 1 else 0
            val synonymous = sampleMutations.filter { x -> x.impact == Impact.SYNONYMOUS }.count()

            val filteredAndSorted = sampleMutations.filter { x -> x.impact != Impact.UNKNOWN && x.impact != Impact.SYNONYMOUS }.sortedWith(oncoComparator)
            if (filteredAndSorted.isEmpty()) {
                return empty.copy(gene = gene, synonymous = synonymous)
            }

            val worst = filteredAndSorted[0]
            val redundant = filteredAndSorted.size - 1
            val isKnown = worst.isHotspot || (worst.impact == Impact.INFRAME && worst.repeatCount < 8)
            val isUnknown = !isKnown
            val result = empty.copy(gene = gene, synonymous = synonymous, redundant = redundant)
            return result.add(worst.impact, isKnown.toInt(), isUnknown.toInt())
        }

        private fun tsgSampleSummary(gene: String, sampleMutations: List<DndsMutation>): MutationsGene {
            fun Boolean.toInt() = if (this) 1 else 0

            val synonymous = sampleMutations.filter { x -> x.impact == Impact.SYNONYMOUS }.count()
            val filteredAndSorted = sampleMutations.filter { x -> x.impact != Impact.UNKNOWN && x.impact != Impact.SYNONYMOUS }.sortedWith(tsgComparator)
            if (filteredAndSorted.isEmpty()) {
                return empty.copy(gene = gene, synonymous = synonymous)
            }

            val worst = filteredAndSorted[0]
            val isKnown = worst.isHotspot || (worst.isBiallelic && worst.impact != Impact.MISSENSE)
            val isUnknown = !isKnown
            val isMultiHit = !isKnown && filteredAndSorted.size > 1
            val redundant = filteredAndSorted.size - 1 - isMultiHit.toInt()

            var result = empty.copy(gene = gene, synonymous = synonymous, redundant = redundant)
            result = result.add(worst.impact, isKnown.toInt(), isUnknown.toInt())
            if (isMultiHit) {
                val secondWorst = filteredAndSorted[1]
                result.add(secondWorst.impact, 0, 1)
            }

            return result
        }

        private fun summary(mutations: List<DndsMutation>, sampleSummary: (String, List<DndsMutation>) -> MutationsGene): List<MutationsGene> {
            val result = mutableListOf<MutationsGene>()
            for (geneMutations in mutations.sortedPartition { x -> x.gene }) {
                val gene = geneMutations[0].gene
                result.add(geneSummary(gene, geneMutations, sampleSummary))
            }

            return result
        }

        private fun geneSummary(gene: String, geneMutations: List<DndsMutation>, sampleSummary: (String, List<DndsMutation>) -> MutationsGene): MutationsGene {
            var result = empty.copy(gene = gene)

            for (sampleMutations in geneMutations.sortedPartition { x -> x.sample }) {
                val sample = sampleMutations[0].sample
                result = result.combine(sampleSummary(gene, sampleMutations))
            }

            return result
        }

        private fun List<DndsMutation>.sortedPartition(key: (DndsMutation) -> String): List<List<DndsMutation>> {
            val result = mutableListOf<List<DndsMutation>>()
            val sorted = this.sortedBy { x -> key(x) }
            var currentKey = ""
            var keyMutations = mutableListOf<DndsMutation>()
            for (mutation in sorted) {
                if (key(mutation) != currentKey) {
                    if (keyMutations.isNotEmpty()) {
                        result.add(keyMutations)
                    }

                    keyMutations = mutableListOf()
                    currentKey = key(mutation)
                }
                keyMutations.add(mutation)
            }

            if (keyMutations.isNotEmpty()) {
                result.add(keyMutations)
            }

            return result
        }
    }

}