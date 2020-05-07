package com.hartwig.hmftools.gripss

import com.google.common.collect.Sets
import htsjdk.variant.variantcontext.GenotypeBuilder
import htsjdk.variant.variantcontext.VariantContext
import htsjdk.variant.variantcontext.VariantContextBuilder
import htsjdk.variant.vcf.VCFCodec
import htsjdk.variant.vcf.VCFHeader
import htsjdk.variant.vcf.VCFHeaderVersion
import org.junit.Assert.assertFalse
import org.junit.Assert.assertTrue
import org.junit.Test

class StructuralVariantContextTest {

    private val codec: VCFCodec by lazy {
        val codec = VCFCodec()
        val header = VCFHeader(Sets.newHashSet(), Sets.newHashSet("NORMAL", "TUMOR"))
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2)
        return@lazy codec
    }

    @Test
    fun testImpreciseFilter() {
        val victim = createBreakEnd()
        assertFalse(victim.toSv().impreciseFilter())
        assertFalse(victim.setAttribute("IMPRECISE", false).toSv().impreciseFilter())
        assertTrue(victim.setAttribute("IMPRECISE", true).toSv().impreciseFilter())
    }

    @Test
    fun testNormalSupportFilter() {
        assertTrue(createBreakEnd().fragmentSupport(3, 9).toSv().normalSupportFilter(0.03))
        assertFalse(createBreakEnd().fragmentSupport(3, 100).toSv().normalSupportFilter(0.03))
        assertTrue(createBreakEnd().fragmentSupport(3, 100).toSv().normalSupportFilter(0.0299))
    }

    @Test
    fun testQualFilter() {
        val breakEnd = createBreakEnd().qual(200).toSv();
        assertTrue(breakEnd.qualFilter(201, 1000))
        assertFalse(breakEnd.qualFilter(200, 1000))
        assertFalse(breakEnd.qualFilter(199, 1000))

        val breakPoint = createBreakPoint().qual(200).toSv();
        assertTrue(breakPoint.qualFilter(1000, 201))
        assertFalse(breakPoint.qualFilter(1000, 200))
        assertFalse(breakPoint.qualFilter(1000, 199))
    }

    @Test
    fun testPolyGFilters() {
        assertFalse(createVariant("A", "A" + "G".repeat(16) + "[1:123[").toSv().polyGCFilter())
        assertTrue(createVariant("A", "A" + "G".repeat(16) + ".").toSv().polyGCFilter())
        assertTrue(createVariant("A", "A" + "C".repeat(16) + ".").toSv().polyGCFilter())
        assertFalse(createVariant("A", "A" + "G".repeat(15) + ".").toSv().polyGCFilter())
        assertFalse(createVariant("A", "A" + "C".repeat(15) + ".").toSv().polyGCFilter())
        assertFalse(createVariant("A", "A" + "A".repeat(16) + ".").toSv().polyGCFilter())
        assertFalse(createVariant("A", "A" + "T".repeat(16) + ".").toSv().polyGCFilter())
    }

    private fun createBreakEnd(): VariantContext {
        return createVariant("C", ".A")
    }

    private fun createBreakPoint(): VariantContext {
        return createVariant("C", "A[2:1000[")
    }

    private fun createVariant(ref: String, alt: String): VariantContext {
        val line = "1\t1\tid1\t${ref}\t${alt}\t100\tPASS\tinfo\tGT:BVF:VF:REF:REFPAIR\t./.:1:1:1:1\t./.:10:10:1:1"
        return codec.decode(line);
    }

    private fun VariantContext.toSv(): StructuralVariantContext = StructuralVariantContext(this)

    private fun VariantContext.qual(qual: Int): VariantContext {
        val builder = VariantContextBuilder(this)
        builder.log10PError(qual.toDouble() / -10.0)
        return builder.make()
    }

    private fun VariantContext.fragmentSupport(normal: Int, tumor: Int): VariantContext {
        return this.addGenotypeAttribute("BVF", normal, tumor).addGenotypeAttribute("VF", normal, tumor);
    }

    private fun VariantContext.setAttribute(attribute: String, value: Any): VariantContext {
        val builder = VariantContextBuilder(this)
        builder.attribute(attribute, value)
        return builder.make();
    }

    private fun VariantContext.addGenotypeAttribute(attribute: String, normal: Int, tumor: Int): VariantContext {
        val normalBuilder = GenotypeBuilder(this.getGenotype(0))
        normalBuilder.attribute(attribute, normal)

        val tumorBuilder = GenotypeBuilder(this.getGenotype(1))
        tumorBuilder.attribute(attribute, tumor)

        val builder = VariantContextBuilder(this)
        builder.genotypes(listOf(normalBuilder.make(), tumorBuilder.make()))
        return builder.make();
    }

}
