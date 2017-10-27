package com.hartwig.hmftools.strelka;

import static org.junit.Assert.assertEquals;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Streams;
import com.google.common.io.Resources;

import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class MNVMergerTest {
    private static final File VCF_FILE = new File(Resources.getResource("mnvs.vcf").getPath());
    private static final VCFFileReader VCF_FILE_READER = new VCFFileReader(VCF_FILE, false);
    private static final List<VariantContext> VARIANTS = Streams.stream(VCF_FILE_READER).collect(Collectors.toList());

    // MIVO: 170756001: (C->T),  170756002: (G->T)
    @Test
    public void correctlyMerges2Variants() {
        final VariantContext mergedVariant =
                MNVMerger.mergeVariants(Lists.newArrayList(VARIANTS.get(11), VARIANTS.get(12)), Maps.newHashMap());
        assertEquals("CG", mergedVariant.getReference().getBaseString());
        assertEquals(1, mergedVariant.getAlternateAlleles().size());
        assertEquals("TT", mergedVariant.getAlternateAllele(0).getBaseString());
    }

    // MIVO: 170756001: (C->T),  170756003: (A->T)    gap: 170756002 (G)
    @Test
    public void correctlyMerges2VariantsWithGap() {
        final Map<Integer, Character> gaps = Maps.newHashMap();
        gaps.put(170756002, 'G');
        final VariantContext mergedVariant = MNVMerger.mergeVariants(Lists.newArrayList(VARIANTS.get(11), VARIANTS.get(13)), gaps);
        assertEquals("CGA", mergedVariant.getReference().getBaseString());
        assertEquals(1, mergedVariant.getAlternateAlleles().size());
        assertEquals("TGT", mergedVariant.getAlternateAllele(0).getBaseString());
    }

    // MIVO: 170756001: (C->T),  170756002: (G->T),  170756003: (A->T)
    @Test
    public void correctlyMerges3Variants() {
        final VariantContext mergedVariant =
                MNVMerger.mergeVariants(Lists.newArrayList(VARIANTS.get(11), VARIANTS.get(12), VARIANTS.get(13)), Maps.newHashMap());
        assertEquals("CGA", mergedVariant.getReference().getBaseString());
        assertEquals(1, mergedVariant.getAlternateAlleles().size());
        assertEquals("TTT", mergedVariant.getAlternateAllele(0).getBaseString());
    }

    // MIVO: 170756001: (C->T),  170756002: (G->T),  170756004: (T->C)
    @Test
    public void correctlyMerges3VariantsWithGap() {
        final Map<Integer, Character> gaps = Maps.newHashMap();
        gaps.put(170756003, 'A');
        final VariantContext mergedVariant =
                MNVMerger.mergeVariants(Lists.newArrayList(VARIANTS.get(11), VARIANTS.get(12), VARIANTS.get(14)), gaps);
        assertEquals("CGAT", mergedVariant.getReference().getBaseString());
        assertEquals(1, mergedVariant.getAlternateAlleles().size());
        assertEquals("TTAC", mergedVariant.getAlternateAllele(0).getBaseString());
    }
}
