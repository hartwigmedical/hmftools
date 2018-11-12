package com.hartwig.hmftools.common.hotspot;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;

public class HotspotEvidenceVCFTest {

    private static final String CHROM = "1";
    private static final String REF = "A";
    private static final int POS = 100;

    private static final int NORMAL_READS = 10;
    private static final int NORMAL_REF_COUNT = 5;
    private static final int TUMOR_READS = 100;
    private static final int TUMOR_REF_COUNT = 40;
    private static final int QUALITY = 500;

    @Test
    public void testCombine() {
        final HotspotEvidence first = create("G", 0, 20);
        final HotspotEvidence second = create("C", 0, 10);
        final HotspotEvidence third = create("T", 1, 30);

        final List<HotspotEvidence> unsorted = Lists.newArrayList(first, third, second);
        Collections.shuffle(unsorted);

        final VariantContext context = new HotspotEvidenceVCF("NORMAL", "TUMOR").create(unsorted);
        assertEquals(NORMAL_READS, context.getGenotype("NORMAL").getDP());
        assertEquals(NORMAL_REF_COUNT, context.getGenotype("NORMAL").getAD()[0]);
        assertEquals(0, context.getGenotype("NORMAL").getAD()[1]);
        assertEquals(0, context.getGenotype("NORMAL").getAD()[2]);
        assertEquals(1, context.getGenotype("NORMAL").getAD()[3]);

        assertEquals(TUMOR_READS, context.getGenotype("TUMOR").getDP());
        assertEquals(TUMOR_REF_COUNT, context.getGenotype("TUMOR").getAD()[0]);
        assertEquals(20, context.getGenotype("TUMOR").getAD()[1]);
        assertEquals(10, context.getGenotype("TUMOR").getAD()[2]);
        assertEquals(30, context.getGenotype("TUMOR").getAD()[3]);
    }

    private static HotspotEvidence create(@NotNull final String alt, int normalAltCount, int tumorAltCount) {
        return create().alt(alt).tumorAltCount(tumorAltCount).normalAltCount(normalAltCount).build();
    }

    private static ImmutableHotspotEvidence.Builder create() {
        return ImmutableHotspotEvidence.builder()
                .chromosome(CHROM)
                .position(POS)
                .ref(REF)
                .alt(REF)
                .qualityScore(QUALITY)
                .type(HotspotEvidenceType.SNV)
                .normalRefCount(NORMAL_REF_COUNT)
                .normalAltCount(0)
                .normalReads(NORMAL_READS)
                .tumorRefCount(TUMOR_REF_COUNT)
                .tumorAltCount(0)
                .tumorReads(TUMOR_READS);
    }

}
