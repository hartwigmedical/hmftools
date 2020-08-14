package com.hartwig.hmftools.paddle.dnds

class DndsMutationComparator(private val includeBiallelicInSort: Boolean) : Comparator<DndsMutation> {
    override fun compare(o1: DndsMutation, o2: DndsMutation): Int {
        fun Boolean.toInt() = if (this) 1 else 0

        if (o1.isHotspot.xor(o2.isHotspot)) {
            return o2.isHotspot.toInt() - o1.isHotspot.toInt()
        }

        if (includeBiallelicInSort && o1.isBiallelic.xor(o2.isBiallelic)) {
            return o2.isBiallelic.toInt() - o1.isBiallelic.toInt()
        }

        return o1.impact.ordinal - o2.impact.ordinal
    }

}