package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.gripss.VariantContextTestFactory.addGenotypeAttribute
import com.hartwig.hmftools.gripss.VariantContextTestFactory.fragmentSupport
import com.hartwig.hmftools.gripss.VariantContextTestFactory.setAttribute
import com.hartwig.hmftools.gripss.VariantContextTestFactory.splitReads
import com.hartwig.hmftools.gripss.VariantContextTestFactory.toSv
import htsjdk.samtools.util.Interval
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import org.junit.Assert.assertFalse
import org.junit.Assert.assertTrue
import org.junit.Test

class StructuralVariantContextTest {

    @Test
    fun testBreakendAssemblyReadPair() {
        val single = sgl();
        assertTrue(single.toSv().breakendAssemblyReadPairsFilter())
        assertFalse(single.setAttribute("BASRP", 1).toSv().breakendAssemblyReadPairsFilter())

        val shortDel = shortDel();
        assertFalse(shortDel.toSv().breakendAssemblyReadPairsFilter())
        assertFalse(shortDel.setAttribute("BASRP", 1).toSv().breakendAssemblyReadPairsFilter())
    }

    @Test
    fun testImpreciseFilter() {
        val victim = sgl()
        assertFalse(victim.toSv().impreciseFilter())
        assertFalse(victim.setAttribute("IMPRECISE", false).toSv().impreciseFilter())
        assertTrue(victim.setAttribute("IMPRECISE", true).toSv().impreciseFilter())
    }

    @Test
    fun testNormalSupportFilter() {
        assertTrue(sgl().fragmentSupport(3, 9).toSv().normalSupportFilter(0.03))
        assertFalse(sgl().fragmentSupport(3, 100).toSv().normalSupportFilter(0.03))
        assertTrue(sgl().fragmentSupport(3, 100).toSv().normalSupportFilter(0.0299))
    }

    @Test
    fun testQualFilter() {
        val breakEnd = sgl().qual(200).toSv();
        assertTrue(breakEnd.qualFilter(201, 1000))
        assertFalse(breakEnd.qualFilter(200, 1000))
        assertFalse(breakEnd.qualFilter(199, 1000))

        val breakPoint = createBreakPoint().qual(200).toSv();
        assertTrue(breakPoint.qualFilter(1000, 201))
        assertFalse(breakPoint.qualFilter(1000, 200))
        assertFalse(breakPoint.qualFilter(1000, 199))
    }

    @Test
    fun testShortDelInsArtifact() {
        val delWithInsertSameSizeAsLength = shortDel(100, 110, "123456789").toSv()
        assertTrue(delWithInsertSameSizeAsLength.shortDelInsertArtifact())

        val delWithInsertLongerSizeAsLength = shortDel(100, 110, "1234567890").toSv()
        assertFalse(delWithInsertLongerSizeAsLength.shortDelInsertArtifact())

        val delWithInsertShorterSizeThanLength = shortDel(100, 110, "12345678").toSv()
        assertFalse(delWithInsertShorterSizeThanLength.shortDelInsertArtifact())
    }

    @Test
    fun testPolyGFilters() {
        val polyGRegion = Interval("1", 1000, 1010)

        assertFalse(createVariant(999, "A", "A" + "G".repeat(100) + "[1:123[").toSv().polyGCFilter(polyGRegion))
        assertTrue(createVariant(1000, "A", "A" + "G".repeat(100) + "[1:123[").toSv().polyGCFilter(polyGRegion))
        assertTrue(createVariant(1010, "A", "A" + "G".repeat(100) + "[1:123[").toSv().polyGCFilter(polyGRegion))
        assertFalse(createVariant(1011, "A", "A" + "G".repeat(100) + "[1:123[").toSv().polyGCFilter(polyGRegion))

        assertFalse(createVariant(2000, "A", "A" + "G".repeat(100) + "[1:999[").toSv().polyGCFilter(polyGRegion))
        assertTrue(createVariant(2000, "A", "A" + "G".repeat(100) + "[1:1000[").toSv().polyGCFilter(polyGRegion))
        assertTrue(createVariant(2000, "A", "A" + "G".repeat(100) + "[1:1010[").toSv().polyGCFilter(polyGRegion))
        assertFalse(createVariant(2000, "A", "A" + "G".repeat(100) + "[1:1011[").toSv().polyGCFilter(polyGRegion))

        assertTrue(createVariant(100, "A", "A" + "G".repeat(16) + ".").toSv().polyGCFilter(polyGRegion))
        assertTrue(createVariant(100, "A", "A" + "C".repeat(16) + ".").toSv().polyGCFilter(polyGRegion))
        assertFalse(createVariant(100, "A", "A" + "G".repeat(15) + ".").toSv().polyGCFilter(polyGRegion))
        assertFalse(createVariant(100, "A", "A" + "C".repeat(15) + ".").toSv().polyGCFilter(polyGRegion))
        assertFalse(createVariant(100, "A", "A" + "A".repeat(16) + ".").toSv().polyGCFilter(polyGRegion))
        assertFalse(createVariant(100, "A", "A" + "T".repeat(16) + ".").toSv().polyGCFilter(polyGRegion))

        assertFalse(createVariant(999, "A", "A.").toSv().polyGCFilter(polyGRegion))
        assertTrue(createVariant(1000, "A", "A.").toSv().polyGCFilter(polyGRegion))
        assertTrue(createVariant(1010, "A", "A.").toSv().polyGCFilter(polyGRegion))
        assertFalse(createVariant(1011, "A", "A.").toSv().polyGCFilter(polyGRegion))
    }

    @Test
    fun testInexactHomologyLengthFilter() {
        val attribute = "IHOMPOS"
        assertFalse(sgl().setAttribute(attribute, listOf(1, 51)).toSv().inexactHomologyLengthFilter(1))

        val bnd = bnd().setAttribute(attribute, listOf(1, 51)).toSv()
        assertFalse(bnd.inexactHomologyLengthFilter(50))
        assertTrue(bnd.inexactHomologyLengthFilter(49))
        assertFalse(bnd.inexactHomologyLengthShortDelFilter(50))
        assertFalse(bnd.inexactHomologyLengthShortDelFilter(49))

        val tinyDel = shortDel(100, 110).setAttribute(attribute, listOf(1, 51)).toSv()
        assertFalse(tinyDel.inexactHomologyLengthFilter(50))
        assertTrue(tinyDel.inexactHomologyLengthFilter(49))
        assertFalse(tinyDel.inexactHomologyLengthShortDelFilter(50))
        assertFalse(tinyDel.inexactHomologyLengthShortDelFilter(49))

        val shortDel = shortDel(100, 310).setAttribute(attribute, listOf(1, 51)).toSv()
        assertFalse(shortDel.inexactHomologyLengthFilter(50))
        assertTrue(shortDel.inexactHomologyLengthFilter(49))
        assertFalse(shortDel.inexactHomologyLengthShortDelFilter(50))
        assertTrue(shortDel.inexactHomologyLengthShortDelFilter(49))

        val shortDup = shortDup().setAttribute(attribute, listOf(1, 51)).toSv()
        assertFalse(shortDup.inexactHomologyLengthFilter(1))
    }

    @Test
    fun testLongDPSupport() {
        assertFalse(sgl().toSv().discordantPairSupportFilter())
        assertFalse(shortDel().toSv().discordantPairSupportFilter())
        assertTrue(bnd().toSv().discordantPairSupportFilter())
        assertFalse(bnd().addGenotypeAttribute("RP", 1, 0).toSv().discordantPairSupportFilter())
        assertFalse(bnd().addGenotypeAttribute("RP", 0, 1).toSv().discordantPairSupportFilter())
        assertFalse(bnd().addGenotypeAttribute("ASRP", 1, 0).toSv().discordantPairSupportFilter())
        assertFalse(bnd().addGenotypeAttribute("ASRP", 0, 1).toSv().discordantPairSupportFilter())
    }

    @Test
    fun testShortSRSupportFilter() {
        val bnd = bnd().splitReads(1, 0).toSv()
        assertFalse(bnd.shortSplitReadTumorFilter())

        assertFalse(shortDel().splitReads(0, 1).toSv().shortSplitReadTumorFilter())
        assertTrue(shortDel().splitReads(0, 0).toSv().shortSplitReadTumorFilter())
    }

    @Test
    fun testShortSRNormalFilter() {
        val bnd = bnd().splitReads(1, 0).toSv()
        assertFalse(bnd.shortSplitReadNormalFilter())

        assertFalse(shortDel().splitReads(0, 1).toSv().shortSplitReadNormalFilter())
        assertTrue(shortDel().splitReads(1, 1).toSv().shortSplitReadNormalFilter())
    }


    @Test
    fun testVariantCreationFunctions() {
        val shortDup = shortDup().toSv()
        assertFalse(shortDup.isSingle)
        assertTrue(shortDup.isShort)
        assertTrue(shortDup.isShortDup)
        assertFalse(shortDup.isShortDel)
        assertFalse(shortDup.isShortIns)

        val shortDel = shortDel().toSv()
        assertFalse(shortDup.isSingle)
        assertTrue(shortDel.isShort)
        assertFalse(shortDel.isShortDup)
        assertTrue(shortDel.isShortDel)
        assertFalse(shortDel.isShortIns)

        val bnd = bnd().toSv()
        assertFalse(bnd.isSingle)
        assertFalse(bnd.isShort)
        assertFalse(bnd.isShortDup)
        assertFalse(bnd.isShortDel)
        assertFalse(bnd.isShortIns)

        val shortInv = shortInv().toSv()
        assertTrue(shortInv.variantType is Inversion)

        assertTrue(sgl().toSv().isSingle)
        assertFalse(createBreakPoint().toSv().isSingle)
    }

    private fun bnd(): VariantContext = createVariant(80, "A", "ATACTGCTACA[2:100[")

    private fun shortInv(position: Int = 100, otherPosition: Int = 110): VariantContext = createVariant(position, "A", "ATT]1:$otherPosition]")

    private fun shortDel(position: Int = 80, otherPosition: Int = 100, insertSequence:String = "TACTGCTACA"): VariantContext = createVariant(position, "A", "A${insertSequence}[1:${otherPosition}[")

    private fun shortDup(position: Int = 110, otherPosition: Int = 100): VariantContext = createVariant(position, "A", "ATACTGCTACA[1:${otherPosition}[")

    private fun sgl(): VariantContext {
        return createVariant(1, "C", ".A")
    }

    private fun createBreakPoint(): VariantContext {
        return createVariant(1, "C", "A[2:1000[")
    }

    private fun createVariant(pos: Int, ref: String, alt: String): VariantContext {
        val line = "1\t${pos}\tid1\t${ref}\t${alt}\t100\tPASS\tinfo\tGT:BVF:VF:REF:REFPAIR\t./.:1:1:1:1\t./.:10:10:1:1"
        return VariantContextTestFactory.decode(line);
    }

    private fun VariantContext.qual(qual: Int): VariantContext {
        val builder = VariantContextBuilder(this)
        builder.log10PError(qual.toDouble() / -10.0)
        return builder.make()
    }


}
