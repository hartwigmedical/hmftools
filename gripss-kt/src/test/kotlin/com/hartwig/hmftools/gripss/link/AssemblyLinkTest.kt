package com.hartwig.hmftools.gripss.link

import com.hartwig.hmftools.gripsskt.StructuralVariantContext
import com.hartwig.hmftools.gripss.VariantContextTestFactory
import com.hartwig.hmftools.gripsskt.store.LinkStore
import com.hartwig.hmftools.gripsskt.link.AssemblyLink
import htsjdk.variant.variantcontext.VariantContext
import junit.framework.Assert.assertEquals
import org.apache.logging.log4j.util.Strings
import org.junit.Test

class AssemblyLinkTest {

    @Test
    fun testNonMatesAreLinked() {
        val variant1 = createVariant(1000, "id1", "T", "]2:1190]T", "mate1", "asm1", "325").toSv()
        val variant2 = createVariant(1010, "id2", "T", "T[3:1190[", "mate2", "asm1", "325").toSv()
        val localLink = LinkStore(AssemblyLink().links(listOf(variant1, variant2)))
        assertEquals("asm1/325", localLink.get("id1"))
        assertEquals("asm1/325", localLink.get("id2"))
    }

    @Test
    fun testNonMatesAreWithDifferentLocationsAreNotLinked() {
        val variant1 = createVariant(1000, "id1", "T", "]2:1190]T", "mate1", "asm1", "325").toSv()
        val variant2 = createVariant(1010, "id2", "T", "T[3:1190[", "mate2", "asm1", "0").toSv()
        val localLink = LinkStore(AssemblyLink().links(listOf(variant1, variant2)))
        assertEquals(Strings.EMPTY, localLink.get("id1"))
        assertEquals(Strings.EMPTY, localLink.get("id2"))
    }

    @Test
    fun testSimpleDelIsNotLinked() {
        val variant1 = createVariant(1000, "id1", "T", "T[1:1010[", "id2", "asm1", "325").toSv()
        val variant2 = createVariant(1010, "id2", "T", "]1:1000]T", "id1", "asm1", "325").toSv()
        val localLink = LinkStore(AssemblyLink().links(listOf(variant1, variant2)))
        assertEquals(Strings.EMPTY, localLink.get("id1"))
        assertEquals(Strings.EMPTY, localLink.get("id2"))
    }

    @Test
    fun testSimpleDelWithSurroundingVariants() {
        val leftOfDel = createVariant(998, "id0", "T", "]1:900]T", "mate0", "asm1", "325").toSv()
        val variant1 = createVariant(1000, "id1", "T", "T[1:1010[", "id2", "asm1", "325").toSv()
        val variant2 = createVariant(1010, "id2", "T", "]1:1000]T", "id1", "asm1", "325").toSv()
        val rightOfDel = createVariant(1012, "id3", "T", "T[1:2010[T", "mate3", "asm1", "325").toSv()
        val localLink = LinkStore(AssemblyLink().links(listOf(leftOfDel, variant1, variant2, rightOfDel)))
        assertEquals("asm1/325-1", localLink.get("id0"))
        assertEquals("asm1/325-1", localLink.get("id1"))
        assertEquals("asm1/325-2", localLink.get("id2"))
        assertEquals("asm1/325-2", localLink.get("id3"))
    }

    fun createVariant(pos: Int, vcfId: String, ref: String, alt: String, mate: String, beid: String, beidl: String): VariantContext {
        val line = "1\t${pos}\t${vcfId}\t${ref}\t${alt}\t100\tPASS\tMATEID=${mate};BEID=${beid};BEIDL=${beidl};AS=10\tGT:BVF:VF:REF:REFPAIR\t./.:1:1:1:1\t./.:10:10:1:1"
        return VariantContextTestFactory.decode(line)
    }

    private fun VariantContext.toSv(): StructuralVariantContext = StructuralVariantContext(this)
}