package com.hartwig.hmftools.ckb.classification;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.CkbTestFactory;

import org.junit.Test;

public class CkbEventAndGeneExtractorTest {

    @Test
    public void canExtractGene() {
        CkbEventAndGeneExtractor extractor = new CkbEventAndGeneExtractor();
        assertEquals("gene", extractor.extractGene(CkbTestFactory.createVariant("gene", "fullName", "variant", null)));
        assertEquals("NTRK2", extractor.extractGene(CkbTestFactory.createVariant("NTRK2", "NTRK2 fusion", "fusion", "fusion")));

        // We fall back on synonyms in case the main symbol that does not exist.
        assertEquals("NTRK2",
                extractor.extractGene(CkbTestFactory.createVariant("does-not-exist",
                        Lists.newArrayList("NTRK2"),
                        "NTRK2 fusion",
                        "fusion",
                        "fusion")));

        // We ignore NO_GENE case
        assertEquals(CkbConstants.NO_GENE,
                extractor.extractGene(CkbTestFactory.createVariant(CkbConstants.NO_GENE, "full_name", "variant", null)));

        // We ignore unmappable genes case
        String firstUnmappableGene = CkbConstants.UNMAPPABLE_GENES.iterator().next();
        assertEquals(firstUnmappableGene,
                extractor.extractGene(CkbTestFactory.createVariant(firstUnmappableGene, "full_name", "variant", null)));
    }

    @Test
    public void canExtractEvent() {
        CkbEventAndGeneExtractor extractor = new CkbEventAndGeneExtractor();
        assertEquals("wild-type", extractor.extractEvent(CkbTestFactory.createVariant("BRAF wild-type", "wild-type", "none")));
        assertEquals("exon 9", extractor.extractEvent(CkbTestFactory.createVariant("KIT exon9", "exon9", "unknown")));
        assertEquals("exon 9", extractor.extractEvent(CkbTestFactory.createVariant("KIT exon 9", "exon 9", "unknown")));
        assertEquals("fusion promiscuous", extractor.extractEvent(CkbTestFactory.createVariant("NTRK2 fusion", "fusion", "fusion")));
        assertEquals("NPM1-ALK", extractor.extractEvent(CkbTestFactory.createVariant("NPM1 - ALK", "NPM1 - ALK", "fusion")));
        assertEquals("unknown", extractor.extractEvent(CkbTestFactory.createVariant("unknown", "unknown", null)));
    }
}