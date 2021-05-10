package com.hartwig.hmftools.ckb.classification;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.CkbTestFactory;

import org.junit.Test;

public class EventAndGeneExtractorTest {

    @Test
    public void canExtractGene() {
        EventAndGeneExtractor extractor = new EventAndGeneExtractor();
        assertEquals("gene", extractor.extractGene(CkbTestFactory.createVariant("gene", "fullName", "variant", null)));
        assertEquals("NPM1 - ALK", extractor.extractGene(CkbTestFactory.createVariant("gene", "NPM1 - ALK", "NPM1 - ALK", "fusion")));
        assertEquals("NTRK2", extractor.extractGene(CkbTestFactory.createVariant("NTRK2", "NTRK2 fusion", "fusion", "fusion")));

        // We fall back on synonyms in case the main symbol that does not exist.
        assertEquals("NTRK2",
                extractor.extractGene(CkbTestFactory.createVariant("does-not-exist",
                        Lists.newArrayList("NTRK2"),
                        "NTRK2 fusion",
                        "fusion",
                        "fusion")));
    }

    @Test
    public void canExtractEvent() {
        EventAndGeneExtractor extractor = new EventAndGeneExtractor();
        assertEquals("wild-type", extractor.extractEvent(CkbTestFactory.createVariant("BRAF wild-type", "wild-type", "none")));
        assertEquals("exon 9", extractor.extractEvent(CkbTestFactory.createVariant("KIT exon9", "exon9", "unknown")));
        assertEquals("exon 9", extractor.extractEvent(CkbTestFactory.createVariant("KIT exon 9", "exon 9", "unknown")));
        assertEquals("fusion promiscuous", extractor.extractEvent(CkbTestFactory.createVariant("NTRK2 fusion", "fusion", "fusion")));
        assertEquals("NPM1-ALK", extractor.extractEvent(CkbTestFactory.createVariant("NPM1 - ALK", "NPM1 - ALK", "fusion")));
        assertEquals("unknown", extractor.extractEvent(CkbTestFactory.createVariant("unknown", "unknown", null)));
    }
}