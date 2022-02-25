package com.hartwig.hmftools.serve.extraction;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.serve.ServeTestFactory;
import com.hartwig.hmftools.serve.actionability.range.ActionableRange;
import com.hartwig.hmftools.serve.actionability.range.ImmutableActionableRange;
import com.hartwig.hmftools.serve.actionability.range.RangeType;
import com.hartwig.hmftools.serve.extraction.codon.KnownCodon;
import com.hartwig.hmftools.serve.extraction.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.extraction.exon.KnownExon;
import com.hartwig.hmftools.serve.extraction.fusion.KnownFusionPair;
import com.hartwig.hmftools.serve.extraction.hotspot.KnownHotspot;
import com.hartwig.hmftools.serve.sources.ImmutableSources;

import org.apache.logging.log4j.util.Strings;
import org.junit.Ignore;
import org.junit.Test;

public class ExtractionFunctionsTest {

    @Test
    @Ignore
    public void canCurateRanges() {
        Set<ActionableRange> actionableRangeSet = Sets.newHashSet();
        ActionableRange range1 = ImmutableActionableRange.builder()
                .from(ServeTestFactory.createTestActionableRange())
                .source(ImmutableSources.builder().source(Knowledgebase.CKB).sourceEvent(Strings.EMPTY).build())
                .gene("BRAF")
                .rank(601)
                .rangeType(RangeType.CODON)
                .start(10)
                .end(20)
                .transcript("transcript")
                .build();
        ActionableRange range2 = ImmutableActionableRange.builder()
                .from(ServeTestFactory.createTestActionableRange())
                .source(ImmutableSources.builder().source(Knowledgebase.CKB).sourceEvent(Strings.EMPTY).build())
                .gene("BRAF")
                .rank(600)
                .rangeType(RangeType.CODON)
                .build();
        actionableRangeSet.add(range1);
        actionableRangeSet.add(range2);

        Set<ActionableRange> actionableRangeSetCorrected = ExtractionFunctions.curationsOfResult(actionableRangeSet);
        List<ActionableRange> listOfActionableRange = new ArrayList<>(actionableRangeSetCorrected);
        ActionableRange actionableRange1 = listOfActionableRange.get(0);
        assertEquals("BRAF", actionableRange1.gene());
        assertEquals(10, actionableRange1.start());
        assertEquals(20, actionableRange1.end());
        assertEquals("transcript", actionableRange1.transcript());

        ActionableRange actionableRange2 = listOfActionableRange.get(1);
        assertEquals("BRAF", actionableRange2.gene());
        assertEquals(140753335, actionableRange2.start());
        assertEquals(140753337, actionableRange2.end());
        assertEquals("ENST00000288602", actionableRange2.transcript());

    }

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