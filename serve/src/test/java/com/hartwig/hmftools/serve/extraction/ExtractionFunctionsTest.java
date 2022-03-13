package com.hartwig.hmftools.serve.extraction;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.ServeTestFactory;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodon;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.extraction.hotspot.KnownHotspot;

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
        KnownHotspot hotspot = merged.knownHotspots().iterator().next();
        assertTrue(hotspot.sources().contains(source1));
        assertTrue(hotspot.sources().contains(source2));

        assertEquals(1, merged.knownCodons().size());
        KnownCodon codon = merged.knownCodons().iterator().next();
        assertTrue(codon.sources().contains(source1));
        assertTrue(codon.sources().contains(source2));

        assertEquals(1, merged.knownExons().size());
        KnownExon exon = merged.knownExons().iterator().next();
        assertTrue(exon.sources().contains(source1));
        assertTrue(exon.sources().contains(source2));

        assertEquals(1, merged.knownCopyNumbers().size());
        KnownCopyNumber copyNumber = merged.knownCopyNumbers().iterator().next();
        assertTrue(copyNumber.sources().contains(source1));
        assertTrue(copyNumber.sources().contains(source2));

        assertEquals(1, merged.knownFusionPairs().size());
        KnownFusionPair fusionPair = merged.knownFusionPairs().iterator().next();
        assertTrue(fusionPair.sources().contains(source1));
        assertTrue(fusionPair.sources().contains(source2));

        assertEquals(2, merged.actionableHotspots().size());
        assertEquals(2, merged.actionableRanges().size());
        assertEquals(2, merged.actionableGenes().size());
        assertEquals(2, merged.actionableFusions().size());
        assertEquals(2, merged.actionableCharacteristics().size());
        assertEquals(2, merged.actionableHLA().size());
    }
}