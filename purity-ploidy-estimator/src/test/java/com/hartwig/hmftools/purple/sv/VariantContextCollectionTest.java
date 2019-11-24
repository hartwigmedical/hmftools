package com.hartwig.hmftools.purple.sv;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class VariantContextCollectionTest {

    private static final String SAMPLE = "sample";

    private VCFCodec codec;

    @Before
    public void setup() {
        codec = createTestCodec();
    }

    @Test
    public void testAddSetsModifiedFlag() {
        final VariantContextCollection victim = new VariantContextCollectionImpl(Lists.newArrayList("1"));
        assertEquals(0, victim.segmentationVariants().size());

        victim.add(codec.decode("1\t192614842\tgridss14_291648b\tT\tTCTCTACACAAG.\t2076.59\tPASS\tSVTYPE=BND\tGT\t./."));
        assertEquals(1, victim.segmentationVariants().size());
    }

    @NotNull
    private static VCFCodec createTestCodec() {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet(SAMPLE));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
        return codec;
    }
}
