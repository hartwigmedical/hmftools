package com.hartwig.hmftools.gripss.store

data class LocationString(val contig: String, val position: Int) {
    companion object {
        operator fun invoke(location: String): LocationString {
            val lastColonIndex: Int = location.lastIndexOf(":")
            val contig = location.substring(0, lastColonIndex)
            val position = location.substring(lastColonIndex + 1, location.length).toInt()
            return LocationString(contig, position)
        }
    }


}