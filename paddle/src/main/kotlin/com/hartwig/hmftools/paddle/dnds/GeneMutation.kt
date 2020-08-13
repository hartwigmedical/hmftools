package com.hartwig.hmftools.paddle.dnds

import org.apache.logging.log4j.LogManager

data class GeneMutation(
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

        private val empty = GeneMutation("empty", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


        fun oncoGeneMutations(mutations: List<DndsMutation>): List<GeneMutation> {
            val result = mutableListOf<GeneMutation>()
            for (geneMutations in mutations.partition { x -> x.gene }) {
                val gene = geneMutations[0].gene
                result.add(oncoGeneMutations(gene, geneMutations))
            }

            return result
        }

        fun oncoGeneMutations(gene: String, geneMutations: List<DndsMutation>): GeneMutation {
            var result = empty.copy(gene = gene)

            for (sampleMutations in geneMutations.partition { x -> x.sample }) {
                val sample = sampleMutations[0].sample
                result = result.add(oncoGeneMutations(sample, gene, sampleMutations))
            }

            return result
        }

        fun oncoGeneMutations(sample: String, gene: String, sampleMutations: List<DndsMutation>): GeneMutation {
            fun Boolean.toInt() = if (this) 1 else 0
            val synonymous = sampleMutations.filter { x -> x.impact == Impact.SYNONYMOUS }.count()

            val filteredAndSorted = sampleMutations.filter { x -> x.impact != Impact.UNKNOWN && x.impact != Impact.SYNONYMOUS }.sorted()
            if (filteredAndSorted.isEmpty()) {
                return empty.copy(gene = gene, synonymous = synonymous)
            }

            val worst = filteredAndSorted[0]
            val redundant = filteredAndSorted.size - 1
            val isKnown = worst.hotspot || (worst.impact == Impact.INFRAME && worst.repeatCount < 8)
            val isUnknown = !isKnown
            return when (worst.impact) {
                Impact.MISSENSE -> empty.copy(gene = gene, synonymous = synonymous, redundant = redundant, knownMissense = isKnown.toInt(), unknownMissense = isUnknown.toInt())
                Impact.INFRAME -> empty.copy(gene = gene, synonymous = synonymous, redundant = redundant, knownInframe = isKnown.toInt(), unknownInframe = isUnknown.toInt())
                Impact.SPLICE -> empty.copy(gene = gene, synonymous = synonymous, redundant = redundant, knownSplice = isKnown.toInt(), unknownSplice = isUnknown.toInt())
                Impact.NONSENSE -> empty.copy(gene = gene, synonymous = synonymous, redundant = redundant, knownNonsense = isKnown.toInt(), unknownNonsense = isUnknown.toInt())
                Impact.FRAMESHIFT -> empty.copy(gene = gene, synonymous = synonymous, redundant = redundant, knownFrameshift = isKnown.toInt(), unknownFrameshift = isUnknown.toInt())
                else -> throw IllegalStateException("Unexpected impact: ${worst.impact}")
            }
        }


        private fun List<DndsMutation>.partition(key: (DndsMutation) -> String): List<List<DndsMutation>> {
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

    fun add(count: GeneMutation): GeneMutation {
        if (gene != count.gene) {
            throw IllegalArgumentException("Incompatible genes: $gene != ${count.gene}")
        }

        return GeneMutation(gene,
                synonymous + count.synonymous,
                redundant + count.redundant,
                knownMissense + count.knownMissense, unknownMissense + count.unknownMissense,
                knownNonsense + count.knownNonsense, unknownNonsense + count.unknownNonsense,
                knownSplice + count.knownSplice, unknownSplice + count.unknownSplice,
                knownInframe + count.knownInframe, unknownInframe + count.unknownInframe,
                knownFrameshift + count.knownFrameshift, unknownFrameshift + count.unknownFrameshift)
    }


}