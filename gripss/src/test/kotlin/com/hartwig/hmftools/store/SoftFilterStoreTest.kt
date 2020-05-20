package com.hartwig.hmftools.store

import com.hartwig.hmftools.gripss.DEDUP
import com.hartwig.hmftools.gripss.MIN_QUAL
import com.hartwig.hmftools.gripss.MIN_SIZE
import com.hartwig.hmftools.gripss.store.SoftFilterStore
import junit.framework.Assert.assertFalse
import junit.framework.Assert.assertTrue
import org.junit.Test

class SoftFilterStoreTest {

    @Test
    fun testEligibleForRescue() {
        val filters: MutableMap<String, Set<String>> = HashMap()
        filters["id1"] = setOf(MIN_QUAL, MIN_SIZE)
        filters["id2"] = setOf(MIN_SIZE, DEDUP)
        val victim = SoftFilterStore(filters)
        assertTrue(victim.isEligibleForRescue("id1"))
        assertFalse(victim.isEligibleForRescue("id2"))
        assertFalse(victim.isEligibleForRescue("id3"))
    }

    @Test
    fun testIsPassing() {
        val filters: MutableMap<String, Set<String>> = HashMap()
        filters["id1"] = setOf(MIN_QUAL, MIN_SIZE)
        filters["id2"] = setOf(MIN_SIZE, DEDUP)
        val victim = SoftFilterStore(filters)
        assertFalse(victim.isPassing("id1"))
        assertFalse(victim.isPassing("id2"))
        assertTrue(victim.isPassing("id3"))
    }

    @Test
    fun testMateFilters() {
        val filters: MutableMap<String, Set<String>> = HashMap()
        filters["id1"] = setOf(MIN_QUAL, MIN_SIZE)
        filters["id2"] = setOf(MIN_SIZE, DEDUP)
        val victim = SoftFilterStore(filters)
        assertTrue(victim.filters("id1", null).isNotEmpty())
        assertTrue(victim.filters("id2", null).isNotEmpty())
        assertTrue(victim.filters("id3", null).isEmpty())
        assertTrue(victim.filters("id3", "id4").isEmpty())
        assertTrue(victim.filters("id3", "id1").isNotEmpty())
    }

}