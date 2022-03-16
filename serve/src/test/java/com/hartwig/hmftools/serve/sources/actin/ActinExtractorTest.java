package com.hartwig.hmftools.serve.sources.actin;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResourceTestFactory;
import com.hartwig.hmftools.serve.sources.actin.classification.ActinClassificationConfig;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinEntry;
import com.hartwig.hmftools.serve.sources.actin.reader.ActinRule;

import org.apache.commons.compress.utils.Lists;
import org.junit.Test;

public class ActinExtractorTest {

    @Test
    public void canExtractFromACTINEntries() {
        EventClassifierConfig config = ActinClassificationConfig.build();
        ActinExtractor extractor = ActinExtractorFactory.buildActinExtractor(config, RefGenomeResourceTestFactory.buildTestResource37());

        List<ActinEntry> actinEntries = Lists.newArrayList();
        actinEntries.add(ActinTestFactory.createTestEntryWithData(ActinRule.AMPLIFICATION_OF_GENE_X, "KIT", ""));
        actinEntries.add(ActinTestFactory.createTestEntryWithData(ActinRule.MUTATION_IN_GENE_X_OF_TYPE_Y, "BRAF", "V600E"));
        actinEntries.add(ActinTestFactory.createTestEntryWithData(ActinRule.FUSION_IN_GENE_X, "NTRK3", ""));
        actinEntries.add(ActinTestFactory.createTestEntryWithData(ActinRule.MUTATION_IN_GENE_X_OF_TYPE_Y, "BRAF", "V600X"));
        actinEntries.add(ActinTestFactory.createTestEntryWithData(ActinRule.MUTATION_IN_GENE_X_OF_TYPE_Y, "BRAF", "exon 2-4"));
        actinEntries.add(ActinTestFactory.createTestEntryWithData(ActinRule.MSI_SIGNATURE, "", "MSI high"));
        actinEntries.add(ActinTestFactory.createTestEntryWithData(ActinRule.TML_OF_AT_LEAST_X, "", "TML >= 450"));
        actinEntries.add(ActinTestFactory.createTestEntryWithData(ActinRule.SPECIFIC_FUSION_OF_X_TO_Y, "", "EML4-ALK"));

        ExtractionResult result = extractor.extract(actinEntries);
        assertEquals(0, result.knownHotspots().size());
        assertEquals(0, result.knownCopyNumbers().size());
        assertEquals(0, result.knownFusionPairs().size());
        assertEquals(1, result.actionableHotspots().size());
        assertEquals(4, result.actionableRanges().size());
        assertEquals(2, result.actionableGenes().size());
        assertEquals(1, result.actionableFusions().size());
        assertEquals(2, result.actionableCharacteristics().size());
    }
}