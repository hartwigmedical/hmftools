package com.hartwig.hmftools.gripss

import htsjdk.samtools.SAMSequenceDictionary

class ContigComparator private constructor(val contigs: Map<String, Int>) : Comparator<String> {

    companion object {

        val defaultContigs = (1..22).map { it.toString() } + "X" + "Y" + "MT" + "M"

        operator fun invoke(dictionary: SAMSequenceDictionary?): ContigComparator {
            val contigs = dictionary?.sequences?.map { it.sequenceName } ?: defaultContigs + defaultContigs.map { "chr$it" }
            return invoke(contigs)
        }

        operator fun invoke(contigs: List<String>): ContigComparator {
            var index = 0
            val map = contigs.associate { Pair(it, index++) }

            if (map.size != contigs.size) {
                throw IllegalArgumentException("Duplicate contigs")
            }

            return ContigComparator(map)
        }
    }


    override fun compare(contig1: String, contig2: String): Int {
        return (this.contigs[contig1] ?: error("contig $contig1 not defined")) - (this.contigs[contig2]
                ?: error("contig $contig2 not defined"))
    }
}