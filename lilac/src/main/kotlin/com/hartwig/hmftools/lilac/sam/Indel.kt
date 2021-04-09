package com.hartwig.hmftools.lilac.sam

data class Indel(val contig: String, val position: Int, val ref: String, val alt: String) {

    companion object {
        fun fromString(line: String): Indel {
            val (contig, position, ref, alt) = line.split("\t")
            return Indel(contig, position.toInt(), ref, alt)
        }
    }

    val isInsert = ref.length < alt.length
    val isDelete = !isInsert
    val length = alt.length - ref.length
    override fun toString(): String {
        return "$contig:$position $ref>$alt"
    }
}