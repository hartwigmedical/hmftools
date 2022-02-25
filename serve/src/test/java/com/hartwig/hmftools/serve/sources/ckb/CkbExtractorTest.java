package com.hartwig.hmftools.serve.sources.ckb;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.ckb.classification.CkbClassificationConfig;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResourceTestFactory;

import org.apache.commons.compress.utils.Lists;
import org.junit.Test;

public class CkbExtractorTest {

    @Test
    public void canExtractFromCKBEntries() {
        EventClassifierConfig config = CkbClassificationConfig.build();
        CkbExtractor extractor = CkbExtractorFactory.buildCkbExtractor(config, RefGenomeResourceTestFactory.buildTestResource37());

        List<CkbEntry> ckbEntries = Lists.newArrayList();
        ckbEntries.add(CkbTestFactory.createEntry("KIT", "amp", "KIT amp", "sensitive", "Actionable", "AB", "cancer", "A", "DOID:162"));
        ckbEntries.add(CkbTestFactory.createEntry("BRAF",
                "V600E",
                "BRAF V600E",
                "sensitive",
                "Actionable",
                "AB",
                "cancer",
                "A",
                "DOID:162"));
        ckbEntries.add(CkbTestFactory.createEntry("NTRK3",
                "fusion promiscuous",
                "NTRK3 fusion promiscuous",
                "sensitive",
                "Actionable",
                "AB",
                "cancer",
                "A",
                "DOID:162"));
        ckbEntries.add(CkbTestFactory.createEntry("BRAF", "V600", "BRAF V600", "sensitive", "Actionable", "AB", "cancer", "A", "DOID:162"));
        ckbEntries.add(CkbTestFactory.createEntry("BRAF",
                "exon 1 deletion",
                "BRAF exon 1 deletion",
                "sensitive",
                "Actionable",
                "AB",
                "cancer",
                "A",
                "DOID:162"));
        ckbEntries.add(CkbTestFactory.createEntry("-", "MSI high", "MSI high", "sensitive", "Actionable", "AB", "cancer", "A", "DOID:162"));
        ckbEntries.add(CkbTestFactory.createEntry("ALk",
                "EML4-ALK",
                "EML4-ALK Fusion",
                "sensitive",
                "Actionable",
                "AB",
                "cancer",
                "A",
                "DOID:162"));

        ExtractionResult result = extractor.extract(ckbEntries);
        assertEquals(1, result.knownHotspots().size());
        assertEquals(1, result.knownCopyNumbers().size());
        assertEquals(1, result.knownFusionPairs().size());
        assertEquals(1, result.actionableHotspots().size());
        assertEquals(2, result.actionableRanges().size());
        assertEquals(2, result.actionableGenes().size());
        assertEquals(1, result.actionableFusions().size());
        assertEquals(1, result.actionableCharacteristics().size());
    }
}