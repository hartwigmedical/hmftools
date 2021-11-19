package com.hartwig.hmftools.gripss.store

class LocationSeek<T : Location> private constructor(private val elements: List<T>) {

    companion object {
        operator fun <T : Location> invoke(elements: List<T>): LocationSeek<T> {
            return LocationSeek(elements.sortedWith(Comparator { x, y -> x.start.compareTo(y.start) }))
        }
    }

    fun any(filter: (T) -> Boolean): Boolean {

        if (elements.isEmpty()) {
            return false
        }

        var loopIndex = 0
        while (loopIndex < elements.size) {
            if (filter(elements[loopIndex])) {
                return true
            }
            loopIndex++
        }

        return false
    }

}