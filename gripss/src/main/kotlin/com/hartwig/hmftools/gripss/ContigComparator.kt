package com.hartwig.hmftools.gripss

import htsjdk.samtools.SAMSequenceDictionary

class ContigComparator private constructor(val contigs: Map<String, Int>) : Comparator<String> {

    companion object {

        operator fun invoke(dictionary: SAMSequenceDictionary) {
            return invoke(dictionary.sequences!!.map { x -> x.sequenceName })
        }

        operator fun invoke(contigs: List<String>) {
            var index = 0
            val map = contigs.associate { Pair(index++, it) }

            if (map.size != contigs.size) {
                throw IllegalArgumentException("Duplicate contigs")
            }

            return ContigComparator(contigs)
        }
    }

    override fun compare(contig1: String, contig2: String): Int {
        return (this.contigs[contig1] ?: error("contig $contig1 not defined")) - (this.contigs[contig2]
                ?: error("contig $contig2 not defined"))
    }
}