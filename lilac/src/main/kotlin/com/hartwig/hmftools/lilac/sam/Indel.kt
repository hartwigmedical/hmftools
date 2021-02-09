package com.hartwig.hmftools.lilac.sam

data class Indel(val contig: String, val position: Int, val ref: String, val alt: String) {
    val isInsert = ref.length < alt.length
    val isDelete = !isInsert
    val length = alt.length - ref.length
    override fun toString(): String {
        return "$contig:$position $ref>$alt"
    }


}