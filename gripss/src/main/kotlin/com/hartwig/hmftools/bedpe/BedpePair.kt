package com.hartwig.hmftools.bedpe

data class BedpeEntry(val contig: String, val start: Int, val end: Int, val orientation: Byte) {

    fun contains(contig: String, start: Int): Boolean {
        return contig == this.contig && start >= this.start && start <= this.end
    }

}

data class BedpePair(val start: BedpeEntry, val end: BedpeEntry) {
    companion object Factory {


        fun fromString(line: String): BedpePair {
            fun String.toOrientation(): Byte = if(this == "+") 1 else -1

            val (contig1, start1, end1, contig2, start2, end2, _, strand1, strand2) = line.split("\t")
            val start = BedpeEntry(contig1, start1.toInt() + 1, end1.toInt(), strand1.toOrientation())
            val end = BedpeEntry(contig2, start2.toInt() + 1, end2.toInt(), strand2.toOrientation())
            return BedpePair(start, end);
        }

        operator fun <T> List<T>.component6(): T = get(5)
        operator fun <T> List<T>.component7(): T = get(6)
        operator fun <T> List<T>.component8(): T = get(7)
        operator fun <T> List<T>.component9(): T = get(8)
    }


}