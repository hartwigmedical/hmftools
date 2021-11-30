package com.hartwig.hmftools.gripsskt.store

import com.hartwig.hmftools.gripsskt.ContigComparator
import com.hartwig.hmftools.gripsskt.StructuralVariantContext

// NOTE: only known fusion pairs are loaded and check - promiscuous locations alone are not

class HotspotStore(private val store: LocationStore) {
    companion object {
        operator fun invoke(compare: ContigComparator, pairedHotspots: List<Breakpoint>): HotspotStore {
            val locationStore = LocationStore.invoke(compare, listOf(), pairedHotspots, 0)
            return HotspotStore(locationStore)
        }

        operator fun invoke(compare: ContigComparator, promiscuousHotspots: List<Breakend>, pairedHotspots: List<Breakpoint>): HotspotStore {
            val locationStore = LocationStore.invoke(compare, promiscuousHotspots, pairedHotspots, 0)
            return HotspotStore(locationStore)
        }
    }

    fun contains(variant: StructuralVariantContext): Boolean {
        if (variant.isSingle) {
            return false
        }

        return containsPromiscuousLeg(variant) || containsPairedHotspot(variant)
    }

    private fun containsPromiscuousLeg(variant: StructuralVariantContext): Boolean {
        return store.contains(variant.startBreakend, 0) || store.contains(variant.endBreakend!!, 0)
    }

    private fun containsPairedHotspot(variant: StructuralVariantContext): Boolean {
        return store.contains(variant.breakpoint!!, 0, 0)
    }

}