package com.hartwig.hmftools.lilac.seq

enum class HlaSequenceIndelType {INSERT, DEL}

sealed class HlaSequenceIndel(val type: HlaSequenceIndelType, val loci: Int, val length: Int) {
    class  HlaSequenceInsert(loci: Int, length: Int): HlaSequenceIndel(HlaSequenceIndelType.INSERT, loci, length)
    class  HlaSequenceDelete(loci: Int, length: Int): HlaSequenceIndel(HlaSequenceIndelType.DEL, loci, length)

    override fun toString(): String {
        return "HlaSequenceIndel(type=$type, loci=$loci, length=$length)"
    }


}
