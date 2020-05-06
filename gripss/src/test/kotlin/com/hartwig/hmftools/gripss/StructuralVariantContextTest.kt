package com.hartwig.hmftools.gripss

import com.google.common.collect.Sets
import htsjdk.variant.variantcontext.GenotypeBuilder
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
        assertTrue(createVariant().addFragmentSupport(3, 9).normalSupportFilter(0.03))
        assertFalse(createVariant().addFragmentSupport(3, 100).normalSupportFilter(0.03))
        assertTrue(createVariant().addFragmentSupport(3, 100).normalSupportFilter(0.0299))
    }

    fun createVariant(): StructuralVariantContext {
        val line = "1\t1\tid1\tC\tA\t100\tPASS\tinfo\tGT\t./.\t./."
        return StructuralVariantContext(codec.decode(line));
    }

    private fun StructuralVariantContext.addFragmentSupport(normal: Int, tumor: Int): StructuralVariantContext {
        return this.addGenotypeAttribute("BVF", normal, tumor).addGenotypeAttribute("VF", normal, tumor);
    }

    private fun StructuralVariantContext.addGenotypeAttribute(attribute: String, normal: Int, tumor: Int): StructuralVariantContext {
        val normalBuilder = GenotypeBuilder(this.context().getGenotype(0))
        normalBuilder.attribute(attribute, normal)

        val tumorBuilder = GenotypeBuilder(this.context().getGenotype(1))
        tumorBuilder.attribute(attribute, tumor)

        val builder = VariantContextBuilder(this.context())
        builder.genotypes(listOf(normalBuilder.make(), tumorBuilder.make()))
        return StructuralVariantContext(builder.make());
    }

}
