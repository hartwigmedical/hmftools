package com.hartwig.hmftools.gripss.dedup

import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.VariantContextTestFactory
import com.hartwig.hmftools.gripss.VariantContextTestFactory.cipos
import com.hartwig.hmftools.gripss.VariantContextTestFactory.setAttribute
import com.hartwig.hmftools.gripss.VariantContextTestFactory.toSv
import com.hartwig.hmftools.gripss.link.Link
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.SoftFilterStore
import com.hartwig.hmftools.gripss.store.VariantStore
import junit.framework.Assert.assertFalse
import junit.framework.Assert.assertTrue
import org.junit.Test

class DedupSingleTest {

    private val fail: Set<String> = setOf("FAIL")
    private val cipos = 10
    private val emptyLinks = LinkStore(listOf())

    @Test
    fun testSingleDedup() {
        val passingSingle = create("sgl", 100, "A.", 1000)
        val passingPair = create("pair", 101, "A[2:222[", 1000)
        val variantStore = VariantStore(listOf(passingSingle, passingPair))
        val softFilterStore = SoftFilterStore(mapOf())

        val victim = DedupSingle(variantStore, softFilterStore, emptyLinks)
        assertFalse(victim.duplicates.contains("pair"))
        assertTrue(victim.duplicates.contains("sgl"))
    }

    @Test
    fun testKeepSingleOverPairBecauseOfQual() {
        val passingSingle = create("sgl", 100, "A.", 1001)
        val passingPair = create("pair", 101, "A[2:222[", 1000)
        val variantStore = VariantStore(listOf(passingSingle, passingPair))
        val softFilterStore = SoftFilterStore(mapOf())

        val victim = DedupSingle(variantStore, softFilterStore, emptyLinks)
        assertFalse(victim.duplicates.contains("sgl"))
        assertTrue(victim.duplicates.contains("pair"))
    }

    @Test
    fun testKeepPairBecauseOfLink() {
        val passingSingle = create("sgl", 100, "A.", 1001)
        val failingPair = create("pair", 101, "A[2:222[", 1000)
        val variantStore = VariantStore(listOf(passingSingle, failingPair))
        val softFilterStore = SoftFilterStore(mapOf(Pair("pair", fail)))
        val link = Link("link", "pair", "other", 0, 0)

        val victim = DedupSingle(variantStore, softFilterStore, LinkStore(listOf(link)))
        assertTrue(victim.duplicates.contains("sgl"))
        assertFalse(victim.duplicates.contains("pair"))
    }

    @Test
    fun testKeepSingleOverPairBecauseOfPairFailing() {
        val passingSingle = create("sgl", 100, "A.", 1)
        val passingPair = create("pair", 101, "A[2:222[", 1000)
        val variantStore = VariantStore(listOf(passingSingle, passingPair))
        val softFilterStore = SoftFilterStore(mapOf(Pair("pair", fail)))

        val victim = DedupSingle(variantStore, softFilterStore, emptyLinks)
        assertFalse(victim.duplicates.contains("sgl"))
        assertTrue(victim.duplicates.contains("pair"))
    }


    @Test
    fun testWithinRange() {
        val passingSingle = create("sgl", 90, "A.", 1000)
        val passingPair = create("pair", 110, "A[2:222[", 1000)
        val variantStore = VariantStore(listOf(passingSingle, passingPair))
        val softFilterStore = SoftFilterStore(mapOf())

        val victim = DedupSingle(variantStore, softFilterStore, emptyLinks)
        assertFalse(victim.duplicates.contains("pair"))
        assertTrue(victim.duplicates.contains("sgl"))
    }

    @Test
    fun testOutsideRange() {
        val passingSingle = create("sgl", 89, "A.", 1000)
        val failingPair = create("pair", 110, "A[2:222[", 1)
        val variantStore = VariantStore(listOf(passingSingle, failingPair))
        val softFilterStore = SoftFilterStore(mapOf(Pair("pair", fail)))

        val victim = DedupSingle(variantStore, softFilterStore, emptyLinks)
        assertFalse(victim.duplicates.contains("pair"))
        assertTrue(victim.duplicates.isEmpty())
    }

    @Test
    fun testMultipleSinglesAndPair() {
        val passingSingle1 = create("sgl1", 100, "A.", 1000)
        val passingSingle2 = create("sgl2", 102, "A.", 1000)
        val failingPair = create("pair", 103, "A[2:222[", 1)
        val variantStore = VariantStore(listOf(passingSingle1, passingSingle2, failingPair))
        val softFilterStore = SoftFilterStore(mapOf(Pair("pair", fail)))

        val victim = DedupSingle(variantStore, softFilterStore, emptyLinks)
        assertFalse(victim.duplicates.contains("pair"))
        assertTrue(victim.duplicates.contains("sgl1"))
        assertTrue(victim.duplicates.contains("sgl2"))
    }

    @Test
    fun testChoosePassingSingle() {
        val passingSingle = create("sgl1", 100, "A.", 1000)
        val failingSingle = create("sgl2", 102, "A.", 1000)
        val variantStore = VariantStore(listOf(passingSingle, failingSingle))
        val softFilterStore = SoftFilterStore(mapOf(Pair("sgl2", fail)))

        val victim = DedupSingle(variantStore, softFilterStore, emptyLinks)
        assertFalse(victim.duplicates.contains("sgl1"))
        assertTrue(victim.duplicates.contains("sgl2"))
    }

    @Test
    fun testChooseHighestQuality() {
        val passingSingle1 = create("sgl1", 100, "A.", 1000)
        val passingSingle2 = create("sgl2", 102, "A.", 1001)
        val variantStore = VariantStore(listOf(passingSingle1, passingSingle2))
        val softFilterStore = SoftFilterStore(mapOf())

        val victim = DedupSingle(variantStore, softFilterStore, emptyLinks)
        assertTrue(victim.duplicates.contains("sgl1"))
        assertFalse(victim.duplicates.contains("sgl2"))
    }

    @Test
    fun testDifferentOrientation() {
        val passingSingle = create("sgl", 100, "A.", 1000)
        val failingPair = create("pair", 100, "[2:222[A", 1)
        val variantStore = VariantStore(listOf(passingSingle, failingPair))
        val softFilterStore = SoftFilterStore(mapOf(Pair("pair", fail)))

        val victim = DedupSingle(variantStore, softFilterStore, emptyLinks)
        assertFalse(victim.duplicates.contains("pair"))
        assertTrue(victim.duplicates.isEmpty())
    }

    @Test
    fun testDoNotUseConfidenceIntervalOfImprecise() {
        val passingSingle = create("sgl", 90, "A.", 1000)
        val imprecise = create("pair", 101, "A[2:222[", 1, precise = false)
        val variantStore = VariantStore(listOf(passingSingle, imprecise))
        val softFilterStore = SoftFilterStore(mapOf(Pair("pair", fail)))

        val victim = DedupSingle(variantStore, softFilterStore, emptyLinks)
        assertFalse(victim.duplicates.contains("pair"))
        assertFalse(victim.duplicates.contains("sgl"))
    }

    @Test
    fun testCanStillMatchAgainstImprecise() {
        val passingSingle = create("sgl", 90, "A.", 1000)
        val imprecise = create("pair", 100, "A[2:222[", 1001, precise = false)
        val variantStore = VariantStore(listOf(passingSingle, imprecise))
        val softFilterStore = SoftFilterStore(mapOf())

        val victim = DedupSingle(variantStore, softFilterStore, emptyLinks)
        assertFalse(victim.duplicates.contains("pair"))
        assertTrue(victim.duplicates.contains("sgl"))
    }

    private fun create(id: String, pos: Int, alt: String, qual: Int, precise: Boolean = true): StructuralVariantContext {
        return VariantContextTestFactory.createVariant("1", pos, id, "A", alt, qual, setOf("."))
                .setAttribute("IMPRECISE", !precise)
                .cipos(Pair(-cipos, cipos), Pair(-cipos, cipos))

                .toSv()
    }
}