package com.hartwig.hmftools.gripss

import htsjdk.samtools.SAMSequenceDictionary

class ContigComparator private constructor(val contigs: Map<String, Int>) : Comparator<String> {

    companion object {

        operator fun invoke(dictionary: SAMSequenceDictionary): ContigComparator {
            return invoke(dictionary.sequences.map { it.sequenceName } )
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

    fun compare(sv1: StructuralVariantContext, sv2: StructuralVariantContext): Int {
        return compare(sv1.contig, sv1.start, sv2.contig, sv2.start)
    }

    fun compare(contig1: String, position1: Int, contig2: String, position2: Int): Int {
        val contigCompare = compare(contig1, contig2)
        return if (contigCompare == 0) {
            position1 - position2
        } else {
            contigCompare
        }
    }

    fun isValidContig(contig: String): Boolean {
        return this.contigs.containsKey(contig)
    }

    override fun compare(contig1: String, contig2: String): Int {
        return (this.contigs[contig1] ?: error("contig $contig1 not defined")) - (this.contigs[contig2]
                ?: error("contig $contig2 not defined"))
    }
}