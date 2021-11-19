package com.hartwig.hmftools.store

import com.hartwig.hmftools.common.gripss.GripssFilters.DEDUP
import com.hartwig.hmftools.common.gripss.GripssFilters.MIN_QUAL
import com.hartwig.hmftools.gripsskt.MIN_LENGTH
import com.hartwig.hmftools.gripsskt.store.SoftFilterStore
import junit.framework.Assert.assertFalse
import junit.framework.Assert.assertTrue
import org.junit.Test

class SoftFilterStoreTest {

    @Test
    fun testIsPassing() {
        val filters: MutableMap<String, Set<String>> = HashMap()
        filters["id1"] = setOf(MIN_QUAL, MIN_LENGTH)
        filters["id2"] = setOf(MIN_LENGTH, DEDUP)
        val victim = SoftFilterStore(filters)
        assertFalse(victim.isPassing("id1"))
        assertFalse(victim.isPassing("id2"))
        assertTrue(victim.isPassing("id3"))
    }

    @Test
    fun testMateFilters() {
        val filters1 = setOf(MIN_QUAL, MIN_LENGTH)
        val filters2 = setOf(MIN_LENGTH, DEDUP)

        val filters: MutableMap<String, Set<String>> = HashMap()
        filters["id1"] = filters1
        filters["id2"] = filters2
        val victim = SoftFilterStore(filters)

        assertTrue(victim.filters("id1", null) == filters1)
        assertTrue(victim.filters("id2", null) == filters2)
        assertTrue(victim.filters("id1", "id2") == filters1 + filters2)
        assertTrue(victim.filters("id3", null).isEmpty())
        assertTrue(victim.filters("id3", "id4").isEmpty())
        assertTrue(victim.filters("id3", "id1") == filters1)
    }

    @Test
    fun containsFilter() {
        val filters: MutableMap<String, Set<String>> = HashMap()
        filters["id1"] = setOf(MIN_QUAL, MIN_LENGTH)
        filters["id2"] = setOf(MIN_LENGTH, DEDUP)
        val victim = SoftFilterStore(filters)

        assertTrue(victim.containsFilter(MIN_QUAL, "id1", null))
        assertTrue(victim.containsFilter(MIN_QUAL, "id1", "id2"))
        assertTrue(victim.containsFilter(MIN_QUAL, "id1", "id3"))

        assertFalse(victim.containsFilter(DEDUP, "id1", null))
        assertTrue(victim.containsFilter(DEDUP, "id1", "id2"))
        assertFalse(victim.containsFilter(DEDUP, "id1", "id3"))
    }

}