package com.hartwig.hmftools.gripss.link

import com.hartwig.hmftools.gripsskt.StructuralVariantContext
import com.hartwig.hmftools.gripss.VariantContextTestFactory
import com.hartwig.hmftools.gripsskt.store.VariantStore
import com.hartwig.hmftools.gripsskt.link.AssemblyLink
import com.hartwig.hmftools.gripsskt.link.DsbLink
import htsjdk.variant.variantcontext.VariantContext
import junit.framework.Assert.assertEquals
import org.junit.Test

class DsbLinkTest {

    @Test
    fun testDsbLink() {
        val variant1 = createVariant(110, "id1", "T", "T[2:100[", "mate1", "asm1", "325").toSv()
        val variant2 = createVariant(120, "id2", "T", "]2:100]T", "mate2", "asm2", "325").toSv()
        val variants = listOf(variant1, variant2)
        val duplicates = setOf<String>()

        val variantStore = VariantStore(variants)
        val assemblyStore = AssemblyLink(variants)

        val dsbLinks = DsbLink(variantStore, assemblyStore, duplicates)
        assertEquals("dsb1", dsbLinks[variant1.vcfId])
        assertEquals("dsb1", dsbLinks[variant2.vcfId])
    }

    @Test
    fun testDsbLinkOneLargeCIPOSAfterInBetweenNonMatching() {
        val variant1 = createVariant(110, "id1", "T", "T[2:100[", "mate1", "asm1", "325").toSv()
        val variant2 = createVariant(220, "id2", "T", "[3:100[T", "mate2", "asm2", "325").toSv()
        val variant3 = createVariant(420, "id2", "T", "]2:100]T", "mate2", "asm2", "325", Pair(-400, 0)).toSv()
        val variants = listOf(variant1, variant2, variant3)
        val duplicates = setOf<String>()

        val variantStore = VariantStore(variants)
        val assemblyStore = AssemblyLink(variants)

        val dsbLinks = DsbLink(variantStore, assemblyStore, duplicates)
        assertEquals("dsb1", dsbLinks[variant1.vcfId])
        assertEquals("dsb1", dsbLinks[variant3.vcfId])
    }

    @Test
    fun testIgnoreAlreadyAssembled() {
        val variant1 = createVariant(110, "id1", "T", "T[2:100[", "mate1", "asm1", "325").toSv()
        val variant2 = createVariant(120, "id2", "T", "]2:100]T", "mate2", "asm1", "325").toSv()
        val variants = listOf(variant1, variant2)
        val duplicates = setOf<String>()

        val variantStore = VariantStore(variants)
        val assemblyStore = AssemblyLink(variants)

        val dsbLinks = DsbLink(variantStore, assemblyStore, duplicates)
        assertEquals("", dsbLinks[variant1.vcfId])
        assertEquals("", dsbLinks[variant2.vcfId])
    }

    @Test
    fun testTooManyCloseByLink() {
        val variant1 = createVariant(110, "id1", "T", ".T", null, "asm1", "325").toSv()
        val variant2 = createVariant(120, "id2", "T", "T.", null, "asm2", "325").toSv()
        val variant3 = createVariant(150, "id3", "T", ".T", null, "asm3", "325").toSv()
        val variants = listOf(variant1, variant2, variant3)
        val duplicates = setOf<String>()

        val variantStore = VariantStore(variants)
        val assemblyStore = AssemblyLink(variants)

        val dsbLinks = DsbLink(variantStore, assemblyStore, duplicates)
        assertEquals("", dsbLinks[variant1.vcfId])
        assertEquals("", dsbLinks[variant2.vcfId])
        assertEquals("", dsbLinks[variant3.vcfId])
    }

    @Test
    fun testIgnoreDuplicates() {
        val variant1 = createVariant(110, "id1", "T", ".T", null, "asm1", "325").toSv()
        val variant2 = createVariant(120, "id2", "T", "T.", null, "asm2", "325").toSv()
        val variant3 = createVariant(150, "id3", "T", ".T", null, "asm3", "325").toSv()
        val variants = listOf(variant1, variant2, variant3)
        val duplicates = setOf(variant1.vcfId)

        val variantStore = VariantStore(variants)
        val assemblyStore = AssemblyLink(variants)

        val dsbLinks = DsbLink(variantStore, assemblyStore, duplicates)
        assertEquals("", dsbLinks[variant1.vcfId])
        assertEquals("dsb1", dsbLinks[variant2.vcfId])
        assertEquals("dsb1", dsbLinks[variant3.vcfId])
    }

    private fun createVariant(pos: Int, vcfId: String, ref: String, alt: String, mate: String?, beid: String, beidl: String): VariantContext {
        val line = "1\t${pos}\t${vcfId}\t${ref}\t${alt}\t100\tPASS\tMATEID=${mate};BEID=${beid};BEIDL=${beidl};AS=10\tGT:BVF:VF:REF:REFPAIR\t./.:1:1:1:1\t./.:10:10:1:1"
        return VariantContextTestFactory.decode(line)
    }

    private fun createVariant(pos: Int, vcfId: String, ref: String, alt: String, mate: String?, beid: String, beidl: String, cipos: Pair<Int, Int>): VariantContext {
        val line = "1\t${pos}\t${vcfId}\t${ref}\t${alt}\t100\tPASS\tMATEID=${mate};BEID=${beid};BEIDL=${beidl};CIPOS=${cipos.first},${cipos.second};AS=10\tGT:BVF:VF:REF:REFPAIR\t./.:1:1:1:1\t./.:10:10:1:1"
        return VariantContextTestFactory.decode(line)
    }

    private fun VariantContext.toSv(): StructuralVariantContext = StructuralVariantContext(this)
}