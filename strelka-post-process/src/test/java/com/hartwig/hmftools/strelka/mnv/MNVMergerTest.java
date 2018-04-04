package com.hartwig.hmftools.strelka.mnv;

import static com.hartwig.hmftools.strelka.StrelkaPostProcessApplication.generateOutputHeader;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

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
import htsjdk.variant.vcf.VCFHeader;

public class MNVMergerTest {
    private static final File VCF_FILE = new File(Resources.getResource("mnvs.vcf").getPath());
    private static final VCFFileReader VCF_FILE_READER = new VCFFileReader(VCF_FILE, false);
    private static final VCFHeader VCF_OUTPUT_HEADER = generateOutputHeader(VCF_FILE_READER.getFileHeader(), "TUMOR");
    private static final MNVMerger MNV_MERGER = ImmutableMNVMerger.of(VCF_OUTPUT_HEADER);
    private static final List<VariantContext> VARIANTS = Streams.stream(VCF_FILE_READER).collect(Collectors.toList());

    // MIVO: 170756001: (C->T),  170756002: (G->T)
    @Test
    public void correctlyMerges2Variants() {
        final VariantContext mergedVariant =
                MNV_MERGER.mergeVariants(Lists.newArrayList(VARIANTS.get(11), VARIANTS.get(12)), Maps.newHashMap());
        assertEquals("CG", mergedVariant.getReference().getBaseString());
        assertEquals(1, mergedVariant.getAlternateAlleles().size());
        assertEquals("TT", mergedVariant.getAlternateAllele(0).getBaseString());
    }

    // MIVO: 170756001: (C->T),  170756003: (A->T)    gap: 170756002 (G)
    @Test
    public void correctlyMerges2VariantsWithGap() {
        final Map<Integer, Character> gaps = Maps.newHashMap();
        gaps.put(170756002, 'G');
        final VariantContext mergedVariant = MNV_MERGER.mergeVariants(Lists.newArrayList(VARIANTS.get(11), VARIANTS.get(13)), gaps);
        assertEquals("CGA", mergedVariant.getReference().getBaseString());
        assertEquals(1, mergedVariant.getAlternateAlleles().size());
        assertEquals("TGT", mergedVariant.getAlternateAllele(0).getBaseString());
    }

    // MIVO: 170756001: (C->T),  170756002: (G->T),  170756003: (A->T)
    @Test
    public void correctlyMerges3Variants() {
        final VariantContext mergedVariant =
                MNV_MERGER.mergeVariants(Lists.newArrayList(VARIANTS.get(11), VARIANTS.get(12), VARIANTS.get(13)), Maps.newHashMap());
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
                MNV_MERGER.mergeVariants(Lists.newArrayList(VARIANTS.get(11), VARIANTS.get(12), VARIANTS.get(14)), gaps);
        assertEquals("CGAT", mergedVariant.getReference().getBaseString());
        assertEquals(1, mergedVariant.getAlternateAlleles().size());
        assertEquals("TTAC", mergedVariant.getAlternateAllele(0).getBaseString());
    }

    // MIVO: 170756001: (C->T),  170756002: (G->T)
    @Test
    public void correctlyMerges2VariantsWithPonCounts() {
        final VariantContext mergedVariant =
                MNV_MERGER.mergeVariants(Lists.newArrayList(VARIANTS.get(18), VARIANTS.get(19)), Maps.newHashMap());
        assertEquals("CG", mergedVariant.getReference().getBaseString());
        assertEquals(1, mergedVariant.getAlternateAlleles().size());
        assertEquals("TT", mergedVariant.getAlternateAllele(0).getBaseString());
        assertEquals(0.5, mergedVariant.getAttribute("MAPPABILITY"));
        assertEquals(10, mergedVariant.getAttribute("GERMLINE_PON_COUNT"));
        assertEquals(7, mergedVariant.getAttribute("SOMATIC_PON_COUNT"));
    }

    // MIVO: 170756001: (C->T),  170756002: (G->T)
    @Test
    public void correctlyMerges2VariantsFirstWithoutAnyPonCount() {
        final VariantContext mergedVariant =
                MNV_MERGER.mergeVariants(Lists.newArrayList(VARIANTS.get(20), VARIANTS.get(21)), Maps.newHashMap());
        assertEquals("CG", mergedVariant.getReference().getBaseString());
        assertEquals(1, mergedVariant.getAlternateAlleles().size());
        assertEquals("TT", mergedVariant.getAlternateAllele(0).getBaseString());
        assertEquals(0.6, mergedVariant.getAttribute("MAPPABILITY"));
        assertNull(mergedVariant.getAttribute("GERMLINE_PON_COUNT"));
        assertNull(mergedVariant.getAttribute("SOMATIC_PON_COUNT"));
    }

    // MIVO: 170756001: (C->T),  170756002: (G->T)
    @Test
    public void correctlyMerges2VariantsSecondWithoutGermlinePonCount() {
        final VariantContext mergedVariant =
                MNV_MERGER.mergeVariants(Lists.newArrayList(VARIANTS.get(22), VARIANTS.get(23)), Maps.newHashMap());
        assertEquals("CG", mergedVariant.getReference().getBaseString());
        assertEquals(1, mergedVariant.getAlternateAlleles().size());
        assertEquals("TT", mergedVariant.getAlternateAllele(0).getBaseString());
        assertEquals(0.7, mergedVariant.getAttribute("MAPPABILITY"));
        assertNull(mergedVariant.getAttribute("GERMLINE_PON_COUNT"));
        assertEquals(5, mergedVariant.getAttribute("SOMATIC_PON_COUNT"));
    }

    // MIVO: 170756001: (C->T),  170756002: (G->T),  170756004: (T->C)
    @Test
    public void correctlyMerges3VariantsWithGapPon() {
        final Map<Integer, Character> gaps = Maps.newHashMap();
        gaps.put(170756003, 'A');
        final VariantContext mergedVariant =
                MNV_MERGER.mergeVariants(Lists.newArrayList(VARIANTS.get(24), VARIANTS.get(25), VARIANTS.get(26)), gaps);
        assertEquals("CGAT", mergedVariant.getReference().getBaseString());
        assertEquals(1, mergedVariant.getAlternateAlleles().size());
        assertEquals("TTAC", mergedVariant.getAlternateAllele(0).getBaseString());
        assertEquals(0.5, mergedVariant.getAttribute("MAPPABILITY"));
        assertNull(mergedVariant.getAttribute("GERMLINE_PON_COUNT"));
        assertNull(mergedVariant.getAttribute("SOMATIC_PON_COUNT"));
    }
}
