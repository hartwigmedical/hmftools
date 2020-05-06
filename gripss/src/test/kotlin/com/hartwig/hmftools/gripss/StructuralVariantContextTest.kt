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
        assertTrue(createSingleBreakEnd().addFragmentSupport(3, 9).toSv().normalSupportFilter(0.03))
        assertFalse(createSingleBreakEnd().addFragmentSupport(3, 100).toSv().normalSupportFilter(0.03))
        assertTrue(createSingleBreakEnd().addFragmentSupport(3, 100).toSv().normalSupportFilter(0.0299))
    }

    private fun createSingleBreakEnd(): VariantContext {
        val line = "1\t1\tid1\tC\tA\t100\tPASS\tinfo\tGT\t./.\t./."
        return codec.decode(line);
    }

    private fun VariantContext.toSv(): StructuralVariantContext = StructuralVariantContext(this)

    private fun VariantContext.addFragmentSupport(normal: Int, tumor: Int): VariantContext {
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
