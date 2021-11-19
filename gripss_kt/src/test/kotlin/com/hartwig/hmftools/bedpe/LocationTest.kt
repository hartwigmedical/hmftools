package com.hartwig.hmftools.bedpe

import com.hartwig.hmftools.gripss.ContigComparator
import com.hartwig.hmftools.gripss.store.Breakend
import com.hartwig.hmftools.gripss.store.Breakpoint
import org.junit.Assert.*
import org.junit.Test

class LocationTest {

    private val defaultContigs = (1..22).map { it.toString() } + "X" + "Y" + "MT" + "M"
    private val contigComparator = ContigComparator(defaultContigs)

    @Test
    fun testDecode() {
        val entry = "1\t9998\t10008\t5\t18606942\t18606952\t.\t9\t-\t-"
        val start = Breakend("1", 9999, 10008, -1)
        val end = Breakend("5", 18606943, 18606952, -1)

        assertEquals(Breakpoint(start, end), Breakpoint.fromBedpe(entry, ContigComparator(defaultContigs)))
    }

    @Test
    fun testBedEquivalence() {
        val victim = Breakend("1", 1000, 1010, -1)
        assertTrue(victim.isEquivalent(victim))
        assertTrue(victim.isEquivalent(victim.copy(start = 1010)))
        assertTrue(victim.isEquivalent(victim.copy(end = 1000)))

        assertFalse(victim.isEquivalent(victim.copy(contig = "2")))
        assertFalse(victim.isEquivalent(victim.copy(orientation = 1)))
        assertFalse(victim.isEquivalent(victim.copy(start = 999, end = 999)))
        assertFalse(victim.isEquivalent(victim.copy(start = 1011, end = 1011)))
    }

    @Test
    fun testReSortOnSameChromosome() {
        val correctEntry = "1\t9997\t9999\t1\t9998\t10008\t.\t9\t-\t+"
        val reverseEntry = "1\t9998\t10008\t1\t9997\t9999\t.\t9\t+\t-"

        assertEquals(Breakpoint.fromBedpe(correctEntry, contigComparator), Breakpoint.fromBedpe(reverseEntry, contigComparator))
    }

    @Test
    fun testReSortOnDifferentChromosome() {
        val correctEntry = "1\t9997\t9999\t3\t9998\t10008\t.\t9\t-\t+"
        val reverseEntry = "3\t9998\t10008\t1\t9997\t9999\t.\t9\t+\t-"

        assertEquals(Breakpoint.fromBedpe(correctEntry, contigComparator), Breakpoint.fromBedpe(reverseEntry, contigComparator))
    }
}