package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.extensions.hasViralSequenceAlignment
import com.hartwig.hmftools.gripss.StructuralVariantContext.Companion.unambiguousNucleotides
import com.hartwig.hmftools.gripss.VariantContextTestFactory.addGenotypeAttribute
import com.hartwig.hmftools.gripss.VariantContextTestFactory.fragmentSupport
import com.hartwig.hmftools.gripss.VariantContextTestFactory.setAttribute
import com.hartwig.hmftools.gripss.VariantContextTestFactory.splitReads
import com.hartwig.hmftools.gripss.VariantContextTestFactory.toSv
import com.hartwig.hmftools.gripss.VariantContextTestFactory.toTumorOnlySv
import htsjdk.samtools.util.Interval
import htsjdk.variant.variantcontext.VariantContext
import org.junit.Assert.*
import org.junit.Test
import kotlin.math.ceil
import kotlin.math.floor

class StructuralVariantContextTest {

    private val defaultContigs = (1..22).map { it.toString() } + "X" + "Y" + "MT" + "M"
    private val contigComparator = ContigComparator(defaultContigs)

    @Test
    fun testUnambiguousNucleotides() {
        val ambiguous = "MGCTABA"
        val fixed = ambiguous.unambiguousNucleotides
        assertEquals("NGCTANA", fixed)
    }

    @Test
    fun testHasViralSequenceAlignment() {
        assertTrue(sgl().setViralSequenceAlignment().hasViralSequenceAlignment(contigComparator))

        assertTrue(sgl().bealn("NC_001526:5225|+|10S415M|60").hasViralSequenceAlignment(contigComparator))
        assertFalse(sgl().bealn("1:5225|+|10S415M|60").hasViralSequenceAlignment(contigComparator))
        assertFalse(sgl().bealn("chr1:5225|+|10S415M|60").hasViralSequenceAlignment(contigComparator))
        assertFalse(sgl().bealn("22:5225|+|10S415M|60").hasViralSequenceAlignment(contigComparator))
        assertFalse(sgl().bealn("chr22:5225|+|10S415M|60").hasViralSequenceAlignment(contigComparator))
        assertTrue(sgl().bealn("23:5225|+|10S415M|60").hasViralSequenceAlignment(contigComparator))
        assertTrue(sgl().bealn("chr23:5225|+|10S415M|60").hasViralSequenceAlignment(contigComparator))
        assertFalse(sgl().bealn("X:5225|+|10S415M|60").hasViralSequenceAlignment(contigComparator))
        assertFalse(sgl().bealn("chrX:5225|+|10S415M|60").hasViralSequenceAlignment(contigComparator))
        assertFalse(sgl().bealn("Y:5225|+|10S415M|60").hasViralSequenceAlignment(contigComparator))
        assertFalse(sgl().bealn("chrY:5225|+|10S415M|60").hasViralSequenceAlignment(contigComparator))
    }

    @Test
    fun testBreakendAssemblyReadPair() {
        val single = sgl().setAttribute("SB", 0.95)
        assertTrue(single.toSv().breakendAssemblyReadPairsFilter())
        assertTrue(single.setAttribute("BASSR", 1).setAttribute("BAQ", 1).toSv().breakendAssemblyReadPairsFilter())
        assertFalse(single.setAttribute("BASRP", 1).toSv().breakendAssemblyReadPairsFilter())
        assertFalse(single.setAttribute("BAQ", 1).toSv().breakendAssemblyReadPairsFilter())

        val shortDel = shortDel()
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
    fun testNormalSupportRelativeFilter() {
        fun assertFilter(expected: Boolean, context: VariantContext, maxNormalRelativeSupport: Double) =
                assertNormalFilter(expected, context) { x -> x.normalSupportRelativeFilter(maxNormalRelativeSupport) }

        assertFilter(true, sgl().fragmentSupport(3, 9), 0.03)
        assertFilter(true, sgl().setViralSequenceAlignment().fragmentSupport(3, 9), 0.03)

        assertFilter(true, bnd().fragmentSupport(3, 9), 0.03)
        assertFilter(true, bnd().setViralSequenceAlignment().fragmentSupport(3, 9), 0.03)

        assertFilter(false, sgl().fragmentSupport(3, 100), 0.03)
        assertFilter(true, sgl().fragmentSupport(3, 100), 0.0299)
    }

    @Test
    fun testNormalSupportAbsoluteFilter() {
        fun assertFilter(expected: Boolean, context: VariantContext, maxNormalAbsoluteSupport: Int) = assertNormalFilter(expected, context) { x -> x.normalSupportAbsoluteFilter(maxNormalAbsoluteSupport) }

        assertFilter(true, sgl().fragmentSupport(4, 0), 3)
        assertFilter(false, sgl().fragmentSupport(3, 0), 3)
        assertFilter(true, sgl().fragmentSupport(3, 0), 2)
    }

    @Test
    fun testTumorQualFilter() {
        val sgl = sgl().qual(200).toSv()
        assertFalse(sgl.tumorQualFilter(200))
        assertTrue(sgl.tumorQualFilter(201))
    }


    @Test
    fun testHardMinQualFilterIsNotRecoveredByHotspot() {
        val config = GripssFilterConfig.default()
        val sgl = sgl().qual(config.hardMinTumorQual - 1).toSv()
        assertTrue(sgl.tumorQualFilter(config.hardMinTumorQual))
        assertFalse(sgl.normalSupportAbsoluteFilter(config.hardMaxNormalAbsoluteSupport))
        assertFalse(sgl.normalSupportRelativeFilter(config.hardMaxNormalRelativeSupport))
        assertTrue(sgl.isHardFilter(config, contigComparator, true))
        assertTrue(sgl.isHardFilter(config, contigComparator, false))
    }

    @Test
    fun testHardMaxNormalAbsoluteSupportIsRecoveredByHotspot() {
        val config = GripssFilterConfig.default()
        val normalSupport = config.hardMaxNormalAbsoluteSupport + 1
        val tumorSupport = ceil(normalSupport / config.softMaxNormalRelativeSupport).toInt() - 1
        val sgl = sgl().qual(config.hardMinTumorQual).fragmentSupport(normalSupport, tumorSupport).toSv()

        assertFalse(sgl.tumorQualFilter(config.hardMinTumorQual))
        assertFalse(sgl.normalSupportRelativeFilter(config.hardMaxNormalRelativeSupport))

        assertTrue(sgl.normalSupportAbsoluteFilter(config.hardMaxNormalAbsoluteSupport))
        assertTrue(sgl.normalSupportRelativeFilter(config.softMaxNormalRelativeSupport))

        assertFalse(sgl.isHardFilter(config, contigComparator, true))
        assertTrue(sgl.isHardFilter(config, contigComparator, false))
    }

    @Test
    fun testHardMaxNormalRequiresSoftRelativeAsWell() {
        val config = GripssFilterConfig.default()
        val normalSupport = config.hardMaxNormalAbsoluteSupport + 1
        val tumorSupport = ceil(normalSupport / config.softMaxNormalRelativeSupport).toInt() + 1
        val sgl = sgl().qual(config.hardMinTumorQual).fragmentSupport(normalSupport, tumorSupport).toSv()
        assertFalse(sgl.tumorQualFilter(config.hardMinTumorQual))
        assertFalse(sgl.normalSupportRelativeFilter(config.hardMaxNormalRelativeSupport))

        assertTrue(sgl.normalSupportAbsoluteFilter(config.hardMaxNormalAbsoluteSupport))
        assertFalse(sgl.normalSupportRelativeFilter(config.softMaxNormalRelativeSupport))

        assertFalse(sgl.isHardFilter(config, contigComparator, true))
        assertFalse(sgl.isHardFilter(config, contigComparator, false))
    }

    @Test
    fun testHardMaxRelativeAbsoluteSupportIsRecoveredByHotspot() {
        val config = GripssFilterConfig.default()
        val normalSupport = config.hardMaxNormalAbsoluteSupport - 1
        val tumorSupport = floor(normalSupport / config.hardMaxNormalRelativeSupport).toInt() - 1
        val sgl = sgl().qual(config.hardMinTumorQual).fragmentSupport(normalSupport, tumorSupport).toSv()
        assertFalse(sgl.tumorQualFilter(config.hardMinTumorQual))
        assertFalse(sgl.normalSupportAbsoluteFilter(config.hardMaxNormalAbsoluteSupport))
        assertTrue(sgl.normalSupportRelativeFilter(config.hardMaxNormalRelativeSupport))

        assertFalse(sgl.isHardFilter(config, contigComparator, true))
        assertTrue(sgl.isHardFilter(config, contigComparator, false))
    }


    @Test
    fun testQualFilter() {
        val breakEnd = sgl().qual(200).toSv()
        assertTrue(breakEnd.qualFilter(201, 1000))
        assertFalse(breakEnd.qualFilter(200, 1000))
        assertFalse(breakEnd.qualFilter(199, 1000))

        val breakPoint = createBreakPoint().qual(200).toSv()
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

        assertFalse(createVariant(999, "A", "A" + "G".repeat(100) + "[1:123[").toSv().polyGCInsertFilter(polyGRegion))
        assertTrue(createVariant(1000, "A", "A" + "G".repeat(100) + "[1:123[").toSv().polyGCInsertFilter(polyGRegion))
        assertTrue(createVariant(1010, "A", "A" + "G".repeat(100) + "[1:123[").toSv().polyGCInsertFilter(polyGRegion))
        assertFalse(createVariant(1011, "A", "A" + "G".repeat(100) + "[1:123[").toSv().polyGCInsertFilter(polyGRegion))

        assertFalse(createVariant(2000, "A", "A" + "G".repeat(100) + "[1:999[").toSv().polyGCInsertFilter(polyGRegion))
        assertTrue(createVariant(2000, "A", "A" + "G".repeat(100) + "[1:1000[").toSv().polyGCInsertFilter(polyGRegion))
        assertTrue(createVariant(2000, "A", "A" + "G".repeat(100) + "[1:1010[").toSv().polyGCInsertFilter(polyGRegion))
        assertFalse(createVariant(2000, "A", "A" + "G".repeat(100) + "[1:1011[").toSv().polyGCInsertFilter(polyGRegion))

        assertTrue(createVariant(100, "A", "A" + "G".repeat(16) + ".").toSv().polyGCInsertFilter(polyGRegion))
        assertTrue(createVariant(100, "A", "A" + "C".repeat(16) + ".").toSv().polyGCInsertFilter(polyGRegion))
        assertFalse(createVariant(100, "A", "A" + "G".repeat(15) + ".").toSv().polyGCInsertFilter(polyGRegion))
        assertFalse(createVariant(100, "A", "A" + "C".repeat(15) + ".").toSv().polyGCInsertFilter(polyGRegion))
        assertFalse(createVariant(100, "A", "A" + "A".repeat(16) + ".").toSv().polyGCInsertFilter(polyGRegion))
        assertFalse(createVariant(100, "A", "A" + "T".repeat(16) + ".").toSv().polyGCInsertFilter(polyGRegion))

        assertFalse(createVariant(999, "A", "A.").toSv().polyGCInsertFilter(polyGRegion))
        assertTrue(createVariant(1000, "A", "A.").toSv().polyGCInsertFilter(polyGRegion))
        assertTrue(createVariant(1010, "A", "A.").toSv().polyGCInsertFilter(polyGRegion))
        assertFalse(createVariant(1011, "A", "A.").toSv().polyGCInsertFilter(polyGRegion))
    }

    @Test
    fun testPolyAFilters() {
        val bnd = bnd()
        assertTrue(bnd.setAttribute("HOMSEQ", "AAAAAAA").toSv().polyATHomologyFilter())
        assertTrue(bnd.setAttribute("HOMSEQ", "TTTTTTT").toSv().polyATHomologyFilter())
        assertFalse(bnd.setAttribute("HOMSEQ", "GTTTTTT").toSv().polyATHomologyFilter())
        assertFalse(bnd.setAttribute("HOMSEQ", "ATTTTTT").toSv().polyATHomologyFilter())
    }

    @Test
    fun testLongDPSupport() {
        fun assertFilter(expected: Boolean, context: VariantContext) = assertNormalFilter(expected, context) { x -> x.discordantPairSupportFilter() }

        assertFilter(false, sgl())
        assertFilter(false, shortDel())
        assertFilter(true, bnd())
        assertFilter(false, bnd().addGenotypeAttribute("RP", 1, 0))
        assertFilter(false, bnd().addGenotypeAttribute("RP", 0, 1))
        assertFilter(false, bnd().addGenotypeAttribute("ASRP", 1, 0))
        assertFilter(false, bnd().addGenotypeAttribute("ASRP", 0, 1))
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
        fun assertFilter(expected: Boolean, context: VariantContext) = assertNormalFilter(expected, context) { x -> x.shortSplitReadNormalFilter() }

        val bnd = bnd().splitReads(1, 0).toSv()
        assertFalse(bnd.shortSplitReadNormalFilter())

        assertFilter(false, shortDel().splitReads(0, 1))
        assertFilter(true, shortDel().splitReads(1, 1))
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

    private fun shortDel(position: Int = 80, otherPosition: Int = 100, insertSequence: String = "TACTGCTACA"): VariantContext = createVariant(position, "A", "A${insertSequence}[1:${otherPosition}[")

    private fun shortDup(position: Int = 110, otherPosition: Int = 100): VariantContext = createVariant(position, "A", "ATACTGCTACA[1:${otherPosition}[")

    private fun sgl(): VariantContext {
        return createVariant(1, "C", ".A")
    }

    private fun createBreakPoint(): VariantContext {
        return createVariant(1, "C", "A[2:1000[")
    }

    private fun createVariant(pos: Int, ref: String, alt: String): VariantContext {
        val line = "1\t${pos}\tid1\t${ref}\t${alt}\t100\tPASS\tinfo\tGT:BVF:VF:REF:REFPAIR\t./.:0:0:1:1\t./.:10:10:1:1"
        val result = VariantContextTestFactory.decode(line)
        return result;
    }

    private fun VariantContext.qual(qual: Int): VariantContext {
        return this.addGenotypeAttribute("QUAL", 0.0, qual.toDouble())
                .addGenotypeAttribute("BQ", 0.0, qual.toDouble())
    }

    private fun VariantContext.setViralSequenceAlignment(): VariantContext {
        return bealn("NC_001526:5225|+|10S415M|60")
    }

    private fun VariantContext.bealn(bealn: String): VariantContext {
        return this.setAttribute("BEALN", bealn)
    }

    private fun assertNormalFilter(expected: Boolean, context: VariantContext, filter: (StructuralVariantContext) -> Boolean) {
        assertEquals(expected, filter(context.toSv()))
        assertFalse(filter(context.toTumorOnlySv()))
    }
}
