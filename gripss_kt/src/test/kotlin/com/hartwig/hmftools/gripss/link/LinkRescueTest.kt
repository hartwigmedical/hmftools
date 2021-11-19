package com.hartwig.hmftools.gripss.link

import com.hartwig.hmftools.gripss.VariantContextTestFactory.createVariant
import com.hartwig.hmftools.gripss.VariantContextTestFactory.toSv
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.SoftFilterStore
import com.hartwig.hmftools.gripss.store.VariantStore
import junit.framework.Assert.assertTrue
import org.junit.Test

class LinkRescueTest {

    private val fail: Set<String> = setOf("FAIL")

    @Test
    fun testPassAcrossLinks() {
        val variantStore = createVariants()

        val link12 = Link("link12", "1start", "2start", 0, 0)
        val link23 = Link("link23", "2end", "3end", 0, 0)
        val link34 = Link("link34", "3start", "4start", 0, 0)
        val linkStore = LinkStore(listOf(link12, link23, link34))

        val softFilterStore = SoftFilterStore(listOf(Pair("1start", fail), Pair("1end", fail), Pair("2start", fail), Pair("2end", fail), Pair("3start", fail), Pair("3end", fail)).toMap())

        val victim = LinkRescue.rescueTransitive(linkStore, softFilterStore, variantStore)
        assertTrue(victim.rescues.contains("1start"))
        assertTrue(victim.rescues.contains("1end"))
        assertTrue(victim.rescues.contains("2start"))
        assertTrue(victim.rescues.contains("2end"))
        assertTrue(victim.rescues.contains("3start"))
        assertTrue(victim.rescues.contains("3end"))
    }

    @Test
    fun testNoInfiniteLoops() {
        val variantStore = createVariants()

        val link12 = Link("link12", "1start", "2start", 0, 0)
        val link23 = Link("link23", "2end", "3end", 0, 0)
        val link34 = Link("link34", "3start", "4start", 0, 0)
        val link41 = Link("link41", "4end", "1end", 0, 0)
        val linkStore = LinkStore(listOf(link12, link23, link34, link41))

        val softFilterStore = SoftFilterStore(listOf(Pair("1start", fail), Pair("1end", fail), Pair("2start", fail), Pair("2end", fail), Pair("3start", fail), Pair("3end", fail)).toMap())

        val victim = LinkRescue.rescueAssembly(linkStore, softFilterStore, variantStore)
        assertTrue(victim.rescues.contains("1start"))
        assertTrue(victim.rescues.contains("1end"))
        assertTrue(victim.rescues.contains("2start"))
        assertTrue(victim.rescues.contains("2end"))
        assertTrue(victim.rescues.contains("3start"))
        assertTrue(victim.rescues.contains("3end"))
    }

    private fun createVariants(): VariantStore {
        val var1 = createVariant("1", 1000, "1start", "A", "A[1:2000[", 100, "1end").toSv()
        val var3 = createVariant("1", 2000, "1end", "A", "]1:1000]A", 100, "1start").toSv()
        val var2 = createVariant("1", 1000, "2start", "A", "A[2:2000[", 100, "2end").toSv()
        val var4 = createVariant("2", 2000, "2end", "A", "]1:1000]A", 100, "2start").toSv()
        val var6 = createVariant("3", 1000, "3start", "A", "A[2:2000[", 100, "3end").toSv()
        val var5 = createVariant("2", 2000, "3end", "A", "]3:1000]A", 100, "3start").toSv()
        val var7 = createVariant("3", 1000, "4start", "A", "A[4:2000[", 100, "4end").toSv()
        val var8 = createVariant("4", 2000, "4end", "A", "]3:1000]A", 100, "4start").toSv()

        return VariantStore(listOf(var1, var2, var3, var4, var5, var6, var7, var8))
    }
}