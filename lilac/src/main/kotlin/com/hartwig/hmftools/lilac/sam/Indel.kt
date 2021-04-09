package com.hartwig.hmftools.lilac.sam

data class Indel(val contig: String, val position: Int, val ref: String, val alt: String) {

    companion object {
        fun fromString(line: String): Indel {
            val (contig, position, ref, alt) = line.split("\t")
            return Indel(contig, position.toInt(), ref, alt)
        }
    }

    fun match(other: Indel): Boolean {
        if (this.contig != other.contig) {
            return false
        }

        if (this.position != other.position) {
            return false
        }

        if (this.ref.length != other.ref.length) {
            return false
        }

        if (this.alt.length != other.alt.length) {
            return false
        }

        if (alt.length > ref.length) {
            // Insert
            if (alt.substring(1) != other.alt.substring(1)) {
                return false
            }
        }

        return true
    }

    val isInsert = ref.length < alt.length
    val isDelete = !isInsert
    val length = alt.length - ref.length
    override fun toString(): String {
        return "$contig:$position $ref>$alt"
    }
}