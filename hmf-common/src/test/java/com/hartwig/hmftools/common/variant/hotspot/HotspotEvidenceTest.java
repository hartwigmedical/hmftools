package com.hartwig.hmftools.common.variant.hotspot;

import static org.junit.Assert.assertEquals;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class HotspotEvidenceTest {

    @Test
    public void testVaf() {
        assertEquals(0.5, create("1", 1, "A", "T").tumorReads(10).tumorAltCount(5).build().vaf(), 0.1);
    }

    @NotNull
    static ImmutableHotspotEvidence.Builder create(@NotNull final String chromosome, long position, String ref, String alt) {
        return ImmutableHotspotEvidence.builder()
                .chromosome(chromosome)
                .position(position)
                .ref(ref)
                .alt(alt)
                .qualityScore(0)
                .normalIndelCount(0)
                .type(HotspotEvidenceType.KNOWN)
                .normalRefCount(0)
                .normalReads(0)
                .tumorRefCount(0)
                .normalAltCount(0)
                .tumorReads(0);
    }
}
