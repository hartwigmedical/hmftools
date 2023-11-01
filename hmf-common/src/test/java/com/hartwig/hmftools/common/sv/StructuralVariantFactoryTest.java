package com.hartwig.hmftools.common.sv;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class StructuralVariantFactoryTest
{
    private static final String SAMPLE = "sample";

    private VCFCodec mCodec;
    private final StructuralVariantFactory mSvFactory = new StructuralVariantFactory(new CompoundFilter(false));

    @Before
    public void setup()
    {
        mCodec = createTestCodec();
    }

    @Test
    public void testSGL()
    {
        final String vcfEntry = "2\t192614842\tgridss14_291648b\tT\tTCTCTACACAAG.\t2076.59\tPASS\tSVTYPE=BND\tGT\t./.";

        final StructuralVariant variant = mSvFactory.createSingleBreakend(mCodec.decode(vcfEntry));
        assertEquals(StructuralVariantType.SGL, variant.type());
    }

    @Test
    public void testUseAFValueOfLeg()
    {
        final String line1 = "2\t321681\tbnd_W\tG\tG]17:198982]\t6\tPASS\tSVTYPE=BND;TAF=1.1,1.2\tGT\t./.";
        final String line2 = "17\t198982\tbnd_W\tG\tG]2:321681]\t6\tPASS\tSVTYPE=BND\tGT\t./.";

        final StructuralVariant variant = mSvFactory.createSV(mCodec.decode(line1), mCodec.decode(line2));
        assertEquals(1.1, variant.start().alleleFrequency(), 0.01);
        assertNull(variant.end().alleleFrequency());
    }

    @Test
    public void testOrientation()
    {
        final String line1 = "2\t321681\tbnd_W\tG\tG]17:198982]\t6\tPASS\tSVTYPE=BND\tGT\t./.";
        final String line2 = "2\t321682\tbnd_V\tT\t]13:123456]T\t6\tPASS\tSVTYPE=BND\tGT\t./.";
        final String line3 = "13\t123456\tbnd_U\tC\tC[2:321682[\t6\tPASS\tSVTYPE=BND\tGT\t./.";
        final String line4 = "13\t123457\tbnd_X\tA\t[17:198983[A\t6\tPASS\tSVTYPE=BND\tGT\t./.";
        final String line5 = "17\t198982\tbnd_Y\tA\tA]2:321681]\t6\tPASS\tSVTYPE=BND\tGT\t./.";
        final String line6 = "17\t198983\tbnd_Z\tC\t[13:123457[C\t6\tPASS\tSVTYPE=BND\tGT\t./.";

        VariantContext c1 = mCodec.decode(line1);
        VariantContext c2 = mCodec.decode(line2);
        VariantContext c3 = mCodec.decode(line3);
        VariantContext c4 = mCodec.decode(line4);
        VariantContext c5 = mCodec.decode(line5);
        VariantContext c6 = mCodec.decode(line6);

        StructuralVariant variant = mSvFactory.createSV(c1, c5);
        final StructuralVariantLeg leg1 = variant.start();
        final StructuralVariantLeg leg5 = variant.end();
        variant = mSvFactory.createSV(c2, c3);
        final StructuralVariantLeg leg2 = variant.start();
        final StructuralVariantLeg leg3 = variant.end();
        variant = mSvFactory.createSV(c4, c6);
        final StructuralVariantLeg leg4 = variant.start();
        final StructuralVariantLeg leg6 = variant.end();

        assertEquals(1, leg1.orientation());
        assertEquals(-1, leg2.orientation());
        assertEquals(1, leg3.orientation());
        assertEquals(-1, leg4.orientation());
        assertEquals(1, leg5.orientation());
        assertEquals(-1, leg6.orientation());
    }

    @NotNull
    private static VCFCodec createTestCodec()
    {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet(SAMPLE));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
        return codec;
    }
}
