package com.hartwig.hmftools.ckb.classification;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.ckb.CkbTestFactory;

import org.junit.Test;

public class CkbEventAndGeneExtractorTest {

    @Test
    public void canExtractGene() {
        assertEquals("gene", CkbEventAndGeneExtractor.extractGene(CkbTestFactory.createVariant("gene", "fullName", "variant", null)));
        assertEquals("NTRK2",
                CkbEventAndGeneExtractor.extractGene(CkbTestFactory.createVariant("NTRK2", "NTRK2 fusion", "fusion", "fusion")));
        assertEquals("EML4 - RET",
                CkbEventAndGeneExtractor.extractGene(CkbTestFactory.createVariant("EML4", "EML4 - RET", "EML4 - RET", "fusion")));

        // We ignore NO_GENE case
        assertEquals(CkbConstants.NO_GENE,
                CkbEventAndGeneExtractor.extractGene(CkbTestFactory.createVariant(CkbConstants.NO_GENE, "full_name", "variant", null)));

        // We ignore unmappable genes case
        String firstUnmappableGene = CkbConstants.UNMAPPABLE_GENES.iterator().next();
        assertEquals(firstUnmappableGene,
                CkbEventAndGeneExtractor.extractGene(CkbTestFactory.createVariant(firstUnmappableGene, "full_name", "variant", null)));
    }

    @Test
    public void canExtractEvent() {
        assertEquals("wild-type",
                CkbEventAndGeneExtractor.extractEvent(CkbTestFactory.createVariant("BRAF wild-type", "wild-type", "none")));
        assertEquals("exon 9", CkbEventAndGeneExtractor.extractEvent(CkbTestFactory.createVariant("KIT exon9", "exon9", "unknown")));
        assertEquals("exon 9", CkbEventAndGeneExtractor.extractEvent(CkbTestFactory.createVariant("KIT exon 9", "exon 9", "unknown")));
        assertEquals("fusion promiscuous",
                CkbEventAndGeneExtractor.extractEvent(CkbTestFactory.createVariant("NTRK2 fusion", "fusion", "fusion")));
        assertEquals("NPM1-ALK", CkbEventAndGeneExtractor.extractEvent(CkbTestFactory.createVariant("NPM1 - ALK", "NPM1 - ALK", "fusion")));
        assertEquals("unknown", CkbEventAndGeneExtractor.extractEvent(CkbTestFactory.createVariant("unknown", "unknown", null)));
    }
}