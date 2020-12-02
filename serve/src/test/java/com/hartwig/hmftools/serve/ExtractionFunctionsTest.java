package com.hartwig.hmftools.serve;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.Knowledgebase;

import org.junit.Test;

public class ExtractionFunctionsTest {

    @Test
    public void canMergeExtractionResults() {
        Knowledgebase source1 = Knowledgebase.VICC_CIVIC;
        Knowledgebase source2 = Knowledgebase.VICC_CGI;
        ExtractionResult result1 = ServeTestFactory.createResultForSource(source1);
        ExtractionResult result2 = ServeTestFactory.createResultForSource(source2);

        ExtractionResult merged = ExtractionFunctions.merge(Lists.newArrayList(result1, result2));

        assertEquals(1, merged.knownHotspots().size());
        assertTrue(merged.knownHotspots().get(0).sources().contains(source1));
        assertTrue(merged.knownHotspots().get(0).sources().contains(source2));

        assertEquals(1, merged.knownCopyNumbers().size());
        assertTrue(merged.knownCopyNumbers().get(0).sources().contains(source1));
        assertTrue(merged.knownCopyNumbers().get(0).sources().contains(source2));

        assertEquals(1, merged.knownFusionPairs().size());
        assertTrue(merged.knownFusionPairs().get(0).sources().contains(source1));
        assertTrue(merged.knownFusionPairs().get(0).sources().contains(source2));

        assertEquals(2, merged.actionableHotspots().size());
        assertEquals(2, merged.actionableRanges().size());
        assertEquals(2, merged.actionableGenes().size());
        assertEquals(2, merged.actionableFusions().size());
        assertEquals(2, merged.actionableSignatures().size());
    }
}