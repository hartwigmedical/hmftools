package com.hartwig.hmftools.paddle.dnds

class DndsMutationComparator(private val knownFunction: (DndsMutation) -> Boolean) : Comparator<DndsMutation> {

    override fun compare(o1: DndsMutation, o2: DndsMutation): Int {
        fun Boolean.toInt() = if (this) 1 else 0

        if (o1.isHotspot.xor(o2.isHotspot)) {
            return o2.isHotspot.toInt() - o1.isHotspot.toInt()
        }

        val o1Known = knownFunction(o1)
        val o2Known = knownFunction(o2)
        if (o1Known.xor(o2Known)) {
            return o2Known.toInt() - o1Known.toInt()
        }

        return o1.impact.ordinal - o2.impact.ordinal
    }
}