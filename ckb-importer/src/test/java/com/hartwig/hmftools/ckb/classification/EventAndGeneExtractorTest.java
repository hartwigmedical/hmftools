package com.hartwig.hmftools.ckb.classification;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.ckb.CkbTestFactory;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class EventAndGeneExtractorTest {

    @Test
    public void canExtractGene() {
        assertGene("gene", CkbTestFactory.createVariant("gene", "fullName", "variant", null));
        assertGene("NPM1 - ALK", CkbTestFactory.createVariant("gene", "NPM1 - ALK", "NPM1 - ALK", "fusion"));
        assertGene("NTRK2", CkbTestFactory.createVariant("NTRK2", "NTRK2 fusion", "fusion", "fusion"));
    }

    @Test
    public void canExtractEvent() {
        assertEvent("wild-type", CkbTestFactory.createVariant("BRAF wild-type", "wild-type", "none"));
        assertEvent("exon 9", CkbTestFactory.createVariant("KIT exon9", "exon9", "unknown"));
        assertEvent("exon 9", CkbTestFactory.createVariant("KIT exon 9", "exon 9", "unknown"));
        assertEvent("fusion promiscuous", CkbTestFactory.createVariant("NTRK2 fusion", "fusion", "fusion"));
        assertEvent("NPM1-ALK", CkbTestFactory.createVariant("NPM1 - ALK", "NPM1 - ALK", "fusion"));
        assertEvent("unknown", CkbTestFactory.createVariant("unknown", "unknown", null));
    }

    private static void assertGene(@NotNull String expected, @NotNull Variant variant) {
        assertEquals(expected, EventAndGeneExtractor.extractGene(variant));
    }

    private static void assertEvent(@NotNull String expected, @NotNull Variant variant) {
        assertEquals(expected, EventAndGeneExtractor.extractEvent(variant));
    }
}