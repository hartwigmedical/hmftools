package com.hartwig.hmftools.paddle.gene

import com.hartwig.hmftools.paddle.dnds.DndsMutation
import com.hartwig.hmftools.paddle.dnds.DndsMutationComparator
import com.hartwig.hmftools.paddle.dnds.Impact
import org.apache.logging.log4j.LogManager

data class GeneMutationSummary(
        val gene: String,
        val synonymous: Int,
        val redundant: Int,
        val knownMissense: Int, val unknownMissense: Int,
        val knownNonsense: Int, val unknownNonsense: Int,
        val knownSplice: Int, val unknownSplice: Int,
        val knownInframe: Int, val unknownInframe: Int,
        val knownFrameshift: Int, val unknownFrameshift: Int) {

    companion object {

        val logger = LogManager.getLogger(this::class.java)

        private val tsgComparator = DndsMutationComparator(true)
        private val oncoComparator = DndsMutationComparator(false)

        private val empty = GeneMutationSummary("empty", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

        fun oncoGeneMutations(mutations: List<DndsMutation>): List<GeneMutationSummary> {
            return summary(mutations, this::oncoSampleSummary)
        }

        fun tsgGeneMutations(mutations: List<DndsMutation>): List<GeneMutationSummary> {
            return summary(mutations, this::tsgSampleSummary)
        }

        private fun oncoSampleSummary(gene: String, sampleMutations: List<DndsMutation>): GeneMutationSummary {
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

        private fun tsgSampleSummary(gene: String, sampleMutations: List<DndsMutation>): GeneMutationSummary {
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

        private fun summary(mutations: List<DndsMutation>, sampleSummary: (String, List<DndsMutation>) -> GeneMutationSummary): List<GeneMutationSummary> {
            val result = mutableListOf<GeneMutationSummary>()
            for (geneMutations in mutations.sortedPartition { x -> x.gene }) {
                val gene = geneMutations[0].gene
                result.add(geneSummary(gene, geneMutations, sampleSummary))
            }

            return result
        }

        private fun geneSummary(gene: String, geneMutations: List<DndsMutation>, sampleSummary: (String, List<DndsMutation>) -> GeneMutationSummary): GeneMutationSummary {
            var result = empty.copy(gene = gene)

            for (sampleMutations in geneMutations.sortedPartition { x -> x.sample }) {
                val sample = sampleMutations[0].sample
                result = result.add(sampleSummary(gene, sampleMutations))
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

    fun add(impact: Impact, known: Int, unknown: Int): GeneMutationSummary {
        return when (impact) {
            Impact.MISSENSE -> copy( synonymous = synonymous, redundant = redundant, knownMissense = known, unknownMissense = unknown)
            Impact.INFRAME -> copy(synonymous = synonymous, redundant = redundant, knownInframe = known, unknownInframe = unknown)
            Impact.SPLICE -> copy( synonymous = synonymous, redundant = redundant, knownSplice = known, unknownSplice = unknown)
            Impact.NONSENSE -> copy(synonymous = synonymous, redundant = redundant, knownNonsense = known, unknownNonsense = unknown)
            Impact.FRAMESHIFT -> copy( synonymous = synonymous, redundant = redundant, knownFrameshift = known, unknownFrameshift = unknown)
            else -> throw IllegalStateException("Unexpected impact: ${impact}")
        }
    }

    fun add(count: GeneMutationSummary): GeneMutationSummary {
        if (gene != count.gene) {
            throw IllegalArgumentException("Incompatible genes: $gene != ${count.gene}")
        }

        return GeneMutationSummary(gene,
                synonymous + count.synonymous,
                redundant + count.redundant,
                knownMissense + count.knownMissense, unknownMissense + count.unknownMissense,
                knownNonsense + count.knownNonsense, unknownNonsense + count.unknownNonsense,
                knownSplice + count.knownSplice, unknownSplice + count.unknownSplice,
                knownInframe + count.knownInframe, unknownInframe + count.unknownInframe,
                knownFrameshift + count.knownFrameshift, unknownFrameshift + count.unknownFrameshift)
    }
}