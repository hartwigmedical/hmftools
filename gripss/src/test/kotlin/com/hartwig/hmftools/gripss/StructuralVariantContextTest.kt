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

    private fun createBreakEnd(): VariantContext {
        val line = "1\t1\tid1\tC\t.A\t100\tPASS\tinfo\tGT:BVF:REF:REFPAIR\t./.:1:1:1\t./.:10:1:1"
        return codec.decode(line);
    }

    private fun createBreakPoint(): VariantContext {
        val line = "1\t1\tid1\tC\tA[2:1000[\t100\tPASS\tinfo\tGT:VF:REF:REFPAIR\t./.:1:1:1\t./.:10:1:1"
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
